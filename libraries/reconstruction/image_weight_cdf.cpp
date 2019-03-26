/***************************************************************************
 * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
 * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "image_weight_cdf.h"
#include "classify_extract_features.h"
#include <data/filters.h>
#include <data/normalize.h>
#include <core/histogram.h>
#include <fstream>

// Read arguments ==========================================================
void ProgWeightCdf::readParams()
{
    fnOdd = getParam("--odd");
    fnEven = getParam("--even");
    fnOut = getParam("-o");

}

// Show ====================================================================
void ProgWeightCdf::showInfo()
{
    if (verbose==0)
        return;
    std::cerr
    << "Odd Image:           " << fnOdd         << std::endl
    << "Even Image:          " << fnEven        << std::endl
    << "Output Image:        " << fnOut       << std::endl;
}

// Usage ===================================================================
void ProgWeightCdf::defineParams()
{
    addUsageLine("Given two aligned micrograph from frames it cast a reliability image");
    addParamsLine("  --odd <micrograph>                       : odd micrograph");
    addParamsLine("  --even <micrograph>                      : even micrograph");
    addParamsLine("  -o <micrograph>                          : Output - reliability image");
}

void ProgWeightCdf::monogenicAmplitude2D(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, FileName fnDebug)
{
	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(freq-freqH);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double iun=DIRECT_MULTIDIM_ELEM(iu,n);
		double un=1.0/iun;
		if (freqH<=un && un<=freq)
		{
			//double H=0.5*(1+cos((un-w1)*ideltal));
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
		} else if (un>freq)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
		}
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	#ifdef DEBUG
	Image<double> filteredvolume;
	filteredvolume = VRiesz;
	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
	#endif



	#ifdef DEBUG
	FileName iternumber;
	iternumber = formatString("_Volume_%i.vol", count);
	Image<double> saveImg2;
	saveImg2() = VRiesz;
	  if (fnDebug.c_str() != "")
	  {
		saveImg2.write(fnDebug+iternumber);
	  }
	saveImg2.clear();
	#endif

	amplitude.resizeNoCopy(VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n) = DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	double uy, ux;
	n=0;
	for(size_t i=0; i<YSIZE(myfftV); ++i)
	{
		uy = VEC_ELEM(freq_fourier,i);
		for(size_t j=0; j<XSIZE(myfftV); ++j)
		{
			ux = VEC_ELEM(freq_fourier,j);
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
			++n;
		}
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	transformer_inv.inverseFourierTransform(fftVRiesz_aux, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
	}

	#ifdef DEBUG
	if (fnDebug.c_str() != "")
	{
	Image<double> saveImg;
	saveImg = amplitude;
	iternumber = formatString("_Amplitude_%i.vol", count);
	saveImg.write(fnDebug+iternumber);
	saveImg.clear();
	}
	#endif // DEBUG
//
	// Low pass filter the monogenic amplitude
	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);
	double raised_w = PI/(freqL-freq);

	n=0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
	{
		double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//		std::cout << "un = " << un << "  freqL = " << freqL << " freq = " << freq << std::endl;
		if ((freqL)>=un && un>=freq)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-freq)));
		}
		else
		{
			if (un>freqL)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
			}
		}
	}
	transformer_inv.inverseFourierTransform();


	#ifdef DEBUG
	saveImg2 = amplitude;
	FileName fnSaveImg2;
	if (fnDebug.c_str() != "")
	{
		iternumber = formatString("_Filtered_Amplitude_%i.vol", count);
		saveImg2.write(fnDebug+iternumber);
	}
	saveImg2.clear();
	#endif // DEBUG
}


/*
void ProgWeightCdf::monogenicAmplitude2D(MultidimArray< std::complex<double> > &myfftV,
		double freq, double freqH, double freqL, MultidimArray<double> &amplitude, int count, FileName fnDebug)
{
	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);
	std::complex<double> J(0,1);

	// Filter the input volume and add it to amplitude
	long n=0;
	double ideltal=PI/(freq-freqH);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
	{
		double iun=DIRECT_MULTIDIM_ELEM(iu,n);
		double un=1.0/iun;
		if (freqH<=un && un<=freq)
		{
			//double H=0.5*(1+cos((un-w1)*ideltal));
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
		} else if (un>freq)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= iun;
		}
	}

	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);

	#ifdef DEBUG
	Image<double> filteredvolume;
	filteredvolume = VRiesz;
	filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
	#endif

	if (fnSpatial!="")
		Vfiltered()=VRiesz;

	#ifdef DEBUG
	FileName iternumber;
	iternumber = formatString("_Volume_%i.vol", count);
	Image<double> saveImg2;
	saveImg2() = VRiesz;
	  if (fnDebug.c_str() != "")
	  {
		saveImg2.write(fnDebug+iternumber);
	  }
	saveImg2.clear();
	#endif

	amplitude.resizeNoCopy(VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate first component of Riesz vector
	double uz, uy, ux;
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				ux = VEC_ELEM(freq_fourier,j);
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second and third components of Riesz vector
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier,i);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}
	transformer_inv.inverseFourierTransform(fftVRiesz, VRiesz);


	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
	{
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
	}
	#ifdef DEBUG
	if (fnDebug.c_str() != "")
	{
	Image<double> saveImg;
	saveImg = amplitude;
	iternumber = formatString("_Amplitude_%i.vol", count);
	saveImg.write(fnDebug+iternumber);
	saveImg.clear();
	}
	#endif // DEBUG
//
	// Low pass filter the monogenic amplitude
	transformer_inv.FourierTransform(amplitude, fftVRiesz, false);
	double raised_w = PI/(freqL-freq);

	n=0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVRiesz)
	{
		double un=1.0/DIRECT_MULTIDIM_ELEM(iu,n);
//		std::cout << "un = " << un << "  freqL = " << freqL << " freq = " << freq << std::endl;
		if ((freqL)>=un && un>=freq)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz,n) *= 0.5*(1 + cos(raised_w*(un-freq)));
		}
		else
		{
			if (un>freqL)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz,n) = 0;
			}
		}
	}
	transformer_inv.inverseFourierTransform();


	#ifdef DEBUG
	saveImg2 = amplitude;
	FileName fnSaveImg2;
	if (fnDebug.c_str() != "")
	{
		iternumber = formatString("_Filtered_Amplitude_%i.vol", count);
		saveImg2.write(fnDebug+iternumber);
	}
	saveImg2.clear();
	#endif // DEBUG
}
*/

void ProgWeightCdf::resolution2eval(int &count_res, double step,
								double &resolution, double &last_resolution,
								double &freq, double &freqL,
								int &last_fourier_idx,
								bool &continueIter,	bool &breakIter,
								bool &doNextIteration)
{
	resolution = maxRes - count_res*step;
	freq = sampling/resolution;
	++count_res;

	double Nyquist = 2*sampling;
	double aux_frequency;
	int fourier_idx;

	DIGFREQ2FFT_IDX(freq, ZSIZE(VRiesz), fourier_idx);

	FFT_IDX2DIGFREQ(fourier_idx, ZSIZE(VRiesz), aux_frequency);

	freq = aux_frequency;

	if (fourier_idx == last_fourier_idx)
	{
		continueIter = true;
		return;
	}

	last_fourier_idx = fourier_idx;
	resolution = sampling/aux_frequency;


	if (count_res == 0)
		last_resolution = resolution;

	if ( ( resolution<Nyquist ))
	{
		breakIter = true;
		return;
	}

	freqL = sampling/(resolution + step);

	int fourier_idx_2;

	DIGFREQ2FFT_IDX(freqL, ZSIZE(VRiesz), fourier_idx_2);

	if (fourier_idx_2 == fourier_idx)
	{
		if (fourier_idx > 0){
			FFT_IDX2DIGFREQ(fourier_idx - 1, ZSIZE(VRiesz), freqL);
		}
		else{
			freqL = sampling/(resolution + step);
		}
	}

}


void ProgWeightCdf::run()
{
	std::cout << "Starting..." << std::endl;
	showInfo();

	Image<double> oddImage, evenImage;

	oddImage.read(fnOdd);
	evenImage.read(fnEven);

	size_t XdimOddImage, XdimEvenImage, YdimOddImage, YdimEvenImage;

	//TODO: make a test to check that images has the same dimensions

	XdimOddImage=XSIZE(oddImage());
	YdimOddImage=YSIZE(oddImage());

	MultidimArray<double> noiseImage, avgImage, realibilityImage, noiseImageforCDF;

	noiseImage.resizeNoCopy(oddImage());
	avgImage.resizeNoCopy(oddImage());

	//Noise model
	noiseImage = oddImage();
	noiseImage -= evenImage();

	//Average Image
	avgImage = oddImage();
	avgImage += evenImage();
	normalize_OldXmipp(noiseImage);
	normalize_OldXmipp(avgImage);



	Image<double> outputResolution;

	outputResolution().initZeros(noiseImage);

	MultidimArray<double> &pOutputResolution = outputResolution();
	MultidimArray<double> amplitudeMS, amplitudeMN;
	MultidimArray<int> &pMask = mask();

	pMask.resizeNoCopy(noiseImage);
	pMask.initConstant(1);


	double significance;
	significance = 0.95;

	double criticalZ=icdf_gauss(significance);
	double criticalW=-1;
	double resolution, resolution_2, last_resolution = 10000;  //A huge value for achieving
												//last_resolution < resolution
	double freq, freqH, freqL;
	double max_meanS = -1e38;
	double cut_value = 0.025;
	int boxsize = 50;

	double freq_step = 1;

	double R_ = freq_step;

	if (R_<0.25)
		R_=0.25;

	sampling = 0.637;

	minRes = 1;
	maxRes = 20;


	double Nyquist = 2*sampling;
	if (minRes<2*sampling)
		minRes = Nyquist;

	bool doNextIteration=true;

	bool lefttrimming = false;
	int last_fourier_idx = -1;

	int count_res = 0;
	FileName fnDebug;

	int iter=0;
	std::vector<double> list;

	std::cout << "Analyzing frequencies" << std::endl;
	std::cout << "                     " << std::endl;
	std::vector<double> noiseValues;
/*
	do
	{
		bool continueIter = false;
		bool breakIter = false;

		resolution2eval(count_res, R_,
						resolution, last_resolution,
						freq, freqH,
						last_fourier_idx, continueIter, breakIter, doNextIteration);

		if (continueIter)
			continue;

		if (breakIter)
			break;

		std::cout << "resolution = " << resolution << std::endl;


		list.push_back(resolution);

		if (iter <2)
			resolution_2 = list[0];
		else
			resolution_2 = list[iter - 2];

		fnDebug = "Signal";

		freqL = freq + 0.01;


		monogenicAmplitude2D(fftV, freq, freqH, freqL, amplitudeMS, iter, fnDebug);
		fnDebug = "Noise";
		monogenicAmplitude2D(*fftN, freq, freqH, freqL, amplitudeMN, iter, fnDebug);


		double sumS=0, sumS2=0, sumN=0, sumN2=0, NN = 0, NS = 0;
		noiseValues.clear();


		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitudeMS)
		{
			if (DIRECT_MULTIDIM_ELEM(pMask, n)>=1)
			{
				double amplitudeValue=DIRECT_MULTIDIM_ELEM(amplitudeMS, n);
				double amplitudeValueN=DIRECT_MULTIDIM_ELEM(amplitudeMN, n);
				sumS  += amplitudeValue;
				sumS2 += amplitudeValue*amplitudeValue;
				noiseValues.push_back(amplitudeValueN);
				sumN  += amplitudeValueN;
				sumN2 += amplitudeValueN*amplitudeValueN;
				++NS;
				++NN;
			}
		}


		#ifdef DEBUG
		std::cout << "NS" << NS << std::endl;
		std::cout << "NVoxelsOriginalMask" << NVoxelsOriginalMask << std::endl;
		std::cout << "NS/NVoxelsOriginalMask = " << NS/NVoxelsOriginalMask << std::endl;
		#endif


			if (NS == 0)
			{
				std::cout << "There are no points to compute inside the mask" << std::endl;
				std::cout << "If the number of computed frequencies is low, perhaps the provided"
						"mask is not enough tight to the volume, in that case please try another mask" << std::endl;
				break;
			}

			double meanS=sumS/NS;
			double sigma2S=sumS2/NS-meanS*meanS;
			double meanN=sumN/NN;
			double sigma2N=sumN2/NN-meanN*meanN;

			// Check local resolution
//			double thresholdNoise;
//			std::sort(noiseValues.begin(),noiseValues.end());
//			thresholdNoise = noiseValues[size_t(noiseValues.size()*significance)];
//			std::cout << "thr Noise " << thresholdNoise << std::endl;


			#ifdef DEBUG
			  std::cout << "Iteration = " << iter << ",   Resolution= " << resolution <<
					  ",   Signal = " << meanS << ",   Noise = " << meanN << ",  Threshold = "
					  << thresholdNoise <<std::endl;
			#endif

			double z=(meanS-meanN)/sqrt(sigma2S/NS+sigma2N/NN);

//			std::cout << "z = " << z << "  zcritical = " << criticalZ << std::endl;

			if (meanS>max_meanS)
				max_meanS = meanS;

			if (meanS<0.001*max_meanS)
			{
				std::cout << "Search of resolutions stopped due to too low signal" << std::endl;
				break;
			}

//			pMask.printShape();
//			amplitudeMS.printShape();
//			pOutputResolution.printShape();


			size_t idx_x_lim, idx_y_lim;

			idx_x_lim = XSIZE(amplitudeMS)/boxsize;
			idx_y_lim = YSIZE(amplitudeMS)/boxsize;


			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(amplitudeMS)
			{
				if (DIRECT_A3D_ELEM(pMask, k,i,j)>=1)
				{
					size_t idx_x = (i/boxsize);
					size_t idx_y = (j/boxsize);

					if (idx_x==idx_x_lim)//9)
						idx_x = idx_x_lim-1;//8;
					if (idx_y==idx_y_lim)//9)
						idx_y = idx_y_lim-1;//8;

					if ( DIRECT_A3D_ELEM(amplitudeMS, k,i,j)>MAT_ELEM(noiseMatrix, idx_x, idx_y) )
					{

						DIRECT_A3D_ELEM(pMask,  k,i,j) = 1;
						DIRECT_A3D_ELEM(pOutputResolution, k,i,j) = resolution;//sampling/freq;
					}
					else{
						DIRECT_A3D_ELEM(pMask,  k,i,j) += 1;
						if (DIRECT_A3D_ELEM(pMask,  k,i,j) >2)
						{
							DIRECT_A3D_ELEM(pMask,  k,i,j) = -1;
							DIRECT_A3D_ELEM(pOutputResolution,  k,i,j) = resolution_2;//maxRes - counter*R_;
						}
					}
				}
			}

			#ifdef DEBUG_MASK
			FileName fnmask_debug;
			fnmask_debug = formatString("maske_%i.vol", iter);
			mask.write(fnmask_debug);
			#endif


			if (doNextIteration)
			{
				if (resolution <= (minRes-0.001))
					doNextIteration = false;
			}

//		}
		iter++;
		last_resolution = resolution;
	} while (doNextIteration);

	Image<double> outputResolutionImage2;
	outputResolutionImage2() = pOutputResolution;
	outputResolutionImage2.write("resultado.vol");

*/

	/*





	realibilityImage=avgImage;

	CDF cdfN;

	noiseImageforCDF = noiseImage;

	cdfN.calculateCDF(noiseImageforCDF);
//	cdfS.calculateCDF(avgImage);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(noiseImage)
	{
		double e=DIRECT_MULTIDIM_ELEM(noiseImage,n);
		double pN=cdfN.getProbability(e);
//			double pS=cdfS.getProbability(e);
//			double pp=pN;
		if (pN>0.05)
		   DIRECT_MULTIDIM_ELEM(realibilityImage,n) = 0; //std::pow(pN,1);//*DIRECT_MULTIDIM_ELEM(avgImage,n);
	}

	Image<double> saveReliableImage;
	saveReliableImage() = realibilityImage;
	saveReliableImage.write(fnOut);

//
//	// Calculate N=Vi-H*S and its energy CDF
//	N().resizeNoCopy(Vi);
//	MultidimArray<double> &mN=N();
//	const MultidimArray<double> &mS=S();
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
//	{
//		double diff=DIRECT_MULTIDIM_ELEM(Vi,n)-DIRECT_MULTIDIM_ELEM(mS,n);
//		DIRECT_MULTIDIM_ELEM(mN,n)=diff*diff;
//	}
//	CDF cdfN;
//	cdfN.calculateCDF(noiseImage);
//
//	// Mask the input volume with a mask value that is proportional to the probability of being larger than noise
//	Vir.resizeNoCopy(Vi);
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mN)
//	{
//		double e=DIRECT_MULTIDIM_ELEM(Vi,n)*DIRECT_MULTIDIM_ELEM(Vi,n);
//		double pN=cdfN.getProbability(e);
//		if (pN<1)
//		{
//			double pS=cdfS.getProbability(e);
//			double pp=pS*pN;
//			DIRECT_MULTIDIM_ELEM(Vir,n)=pp*DIRECT_MULTIDIM_ELEM(Vi,n);
//		}
//		else
//			DIRECT_MULTIDIM_ELEM(Vir,n)=DIRECT_MULTIDIM_ELEM(Vi,n);
//	}
//
//
*/

}
