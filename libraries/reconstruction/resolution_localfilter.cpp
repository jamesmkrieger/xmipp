/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
 * 			   Carlos Oscar S. Sorzano            coss@cnb.csic.es (2018)
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

#include "resolution_localfilter.h"
//#define DEBUG
//#define DEBUG_MASK

void ProgResLocalFilter::readParams()
{
	fnVol = getParam("--vol");
	fnRes = getParam("--resvol");
	fnOut = getParam("-o");
	fnFilt = getParam("--filteredMap");
	sampling = getDoubleParam("--sampling_rate");
	freq_step = getDoubleParam("--step");
	onlyLocalfilter = checkParam("--nosharpening");
	nthrs = getIntParam("--threads");
}


void ProgResLocalFilter::defineParams()
{
	addUsageLine("This function performs a local filter of a map/tomogram based on the local resolution values");
	addParamsLine("  --vol <vol_file=\"\">   			        : Volume");
	addParamsLine("  --resvol <vol_file=\"\">					: Resolution map");
	addParamsLine("  -o <output=\"MGresolution.vol\">			: Local resolution volume (in Angstroms)");
	addParamsLine("  --filteredMap <output=\"filteredMap.vol\">	: Local resolution volume filtered (in Angstroms)");
	addParamsLine("  [--sampling_rate <s=1>]   					: Sampling rate (A/px)");
	addParamsLine("  [--nosharpening]       					: Do not apply sharpnening, only it only applies a local filter"
				  " according to local resolution values");
	addParamsLine("  [--step <s=0.25>]       					: The resolution is computed at a number of frequencies between minimum and");
	addParamsLine("  [--threads <s=4>]               			: Number of threads");
}


void ProgResLocalFilter::produceSideInfo()
{
	std::cout << "Starting..." << std::endl;
	std::cout << "           " << std::endl;

	Image<double> V;

	V.read(fnVol);
	V().setXmippOrigin();

	std::cout << "Map read" << std::endl;
	MultidimArray<double> &inputVol = V();

	//It applies a smoothing of 10 px in the image borders
	int N_smoothing = 10;

	int siz_z = ZSIZE(inputVol)*0.5;
	int siz_y = YSIZE(inputVol)*0.5;
	int siz_x = XSIZE(inputVol)*0.5;

	//Smoothing the boundaries of the volume
	int limit_distance_x = (siz_x-N_smoothing);
	int limit_distance_y = (siz_y-N_smoothing);
	int limit_distance_z = (siz_z-N_smoothing);

	double uz, uy, ux, uz2, u2, uz2y2;
	long n=0;
	n=0;
	for(int k=0; k<ZSIZE(inputVol); ++k)
	{
		uz = (k - siz_z);
		for(int i=0; i<YSIZE(inputVol); ++i)
		{
			uy = (i - siz_y);
			for(int j=0; j<XSIZE(inputVol); ++j)
			{
				ux = (j - siz_x);

				if (abs(ux)>=limit_distance_x)
				{
					DIRECT_MULTIDIM_ELEM(inputVol, n) *= 0.5*(1+cos(PI*(limit_distance_x - abs(ux))/(N_smoothing)));
				}
				if (abs(uy)>=limit_distance_y)
				{
					DIRECT_MULTIDIM_ELEM(inputVol, n) *= 0.5*(1+cos(PI*(limit_distance_y - abs(uy))/(N_smoothing)));
				}
				if (abs(uz)>=limit_distance_z)
				{
					DIRECT_MULTIDIM_ELEM(inputVol, n) *= 0.5*(1+cos(PI*(limit_distance_z - abs(uz))/(N_smoothing)));
				}
				++n;
			}
		}
	}

	FourierTransformer transformer;
	transformer.FourierTransform(inputVol, fftV);

	MultidimArray<double> inputVol2;
	size_t Zdim, Ydim, Xdim, Ndim;
	fftV.getDimensions(Xdim, Ydim, Zdim, Ndim);
//	inputVol2.resizeNoCopy(Ndim, Zdim, Ydim, Xdim);
//
//
//	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftV)
//		DIRECT_MULTIDIM_ELEM(inputVol2, n) = log(fabs(DIRECT_MULTIDIM_ELEM(fftV, n)));
//
//	Image<double> filteredvolume;
//	filteredvolume = inputVol2;
//	filteredvolume.write("fourier.vol");
	iu.resizeNoCopy(Ndim, Zdim, Ydim, Xdim);


	// Calculate frequency map
	n=0;
	for(size_t k=0; k<ZSIZE(fftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
		uz2=uz*uz;

		for(size_t i=0; i<YSIZE(fftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
			uz2y2=uz2+uy*uy;

			for(size_t j=0; j<XSIZE(fftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol),ux);
				u2=uz2y2+ux*ux;
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
				++n;
			}
		}
	}
	V.clear();


	if ( (ZSIZE(fftV) <= XSIZE(fftV)) && (ZSIZE(fftV) <= YSIZE(fftV)) )
		len = ZSIZE(fftV);
	if ( (YSIZE(fftV) <= ZSIZE(fftV)) && (YSIZE(fftV) <= XSIZE(fftV)) )
		len = YSIZE(fftV);
	else
		len = XSIZE(fftV);

	freq_fourier.initZeros(len);


	resVol.read(fnRes);
	std::cout << "resolution read" << std::endl;
	MultidimArray<double> &pResVol = resVol();
	maxRes = -1;
	minRes = 1e38;
	double sumres = 0, sumres2 = 0;

	long counter = 0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(pResVol)
	{
		double resVal = DIRECT_MULTIDIM_ELEM(pResVol, n);

		sumres += resVal;
		sumres2 += resVal*resVal;

		if (resVal>=maxRes)
			maxRes=resVal;
		if (resVal<minRes)
			minRes=resVal;
		++counter;
	}
	double mean = sumres/counter;
	sigma = sumres2/counter-mean*mean;
	std::cout << "sigma = " << sigma << std::endl;
	double u;

	VEC_ELEM(freq_fourier,0) = 1e-38;
	for(size_t k=0; k<len; ++k)
	{
		FFT_IDX2DIGFREQ(k,len, u);
		VEC_ELEM(freq_fourier,k) = u;
	}
}

void ProgResLocalFilter::localfilteredMap(MultidimArray< std::complex<double> > &myfftV,
        MultidimArray<double> &localfilteredVol,
        double &minRes, double &maxRes, double &step)
{

	//Determining the frequency range;
	int lowIdx, highIdx;
	double lowestfreq, highestfreq, freqL, freqH, freq;
	lowestfreq = sampling/maxRes;
	highestfreq = sampling/minRes;

	std::cout << "max = " << lowestfreq << "  " << maxRes << std::endl;
	std::cout << "min = " << highestfreq << "  " << minRes << std::endl;

	DIGFREQ2FFT_IDX(lowestfreq, len, lowIdx);
	DIGFREQ2FFT_IDX(highestfreq, len, highIdx);

	MultidimArray< std::complex<double> > fftVaux;

	MultidimArray<double> filteredMap, filtered_aux;
	filtered_aux.resizeNoCopy(resVol());
	filteredMap.initZeros(resVol());

	std::cout << "freq = " << highIdx << std::endl;
	std::cout << "freq = " << lowIdx << std::endl;

	double lastRes=1e38;


	for (double idx = lowIdx; idx < highIdx; idx++)
	{
		FFT_IDX2DIGFREQ(idx, len, freq);

		if (lastRes - sampling/freq < step)
			continue;

		std::cout << "freq = " << freq << std::endl;
		std::cout << "res = " << sampling/freq << std::endl;

		freqL = freq - 0.02;
		if (freqL < 0)
			freqL = 0.001;

		freqH = freq + 0.02;

		if (freqH > 0.5)
			freqH = 0.5;

		fftVaux.initZeros(fftV);

		// Filter the input volume
		long n=0;
		double ideltal=PI/(freq-freqH);
		double ideltah=PI/(freqL-freq);


		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVaux)
		{
			double iun=DIRECT_MULTIDIM_ELEM(iu,n);
			double un=1.0/iun;
			if (un <= freqH && freq <= un)
			{
				//double H=0.5*(1+cos((un-w1)*ideltal));
				DIRECT_MULTIDIM_ELEM(fftVaux, n) = DIRECT_MULTIDIM_ELEM(fftV, n);
				DIRECT_MULTIDIM_ELEM(fftVaux, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
			}
			if (freqL<=un && un<=freq)
			{
				//double H=0.5*(1+cos((un-w1)*ideltal));
				DIRECT_MULTIDIM_ELEM(fftVaux, n) = DIRECT_MULTIDIM_ELEM(fftV, n);
				DIRECT_MULTIDIM_ELEM(fftVaux, n) *= 0.5*(1+cos((un-freq)*ideltah));//H;
			}
			if (freqL<=un && un<=freqH)
			{
				//double H=0.5*(1+cos((un-w1)*ideltal));
				DIRECT_MULTIDIM_ELEM(fftVaux, n) = DIRECT_MULTIDIM_ELEM(fftV, n);
			}
		}

		transformer_inv.inverseFourierTransform(fftVaux, filtered_aux);

		double stdres = sampling/sqrt(sigma);

		//
		double sumfreq = 0, sumfreq2 = 0;
		long counter = 0;
		//

		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filtered_aux)
		{
			double weight, sumweight, digitalfreq;
			//TODO: make the converstion in the produdesideinfo
			digitalfreq = sampling/DIRECT_MULTIDIM_ELEM(resVol, n);
			counter = 0;
			for (double idx = lowIdx; idx < highIdx; idx++)
				{
				double freq_aux;
				FFT_IDX2DIGFREQ(idx, len, freq_aux);
				sumfreq += freq_aux;
				sumfreq2 += freq_aux*freq_aux;
				}

			double mean = sumfreq/counter;
			stdres = sumfreq2/counter-mean*mean;

			weight = exp(-(digitalfreq - freq)*(digitalfreq - freq)/(stdres));
			DIRECT_MULTIDIM_ELEM(filteredMap, n) += weight*DIRECT_MULTIDIM_ELEM(filtered_aux, n);
			sumweight += weight;
		}
		lastRes = sampling/freq;
	}
}




void ProgResLocalFilter::run()
{
	produceSideInfo();

	//Determining the frequency range;
	int lowIdx, highIdx;
	double lowestfreq, highestfreq, freqL, freqH, freq;
	lowestfreq = sampling/maxRes;
	highestfreq = sampling/minRes;

	std::cout << "max = " << lowestfreq << "  " << maxRes << std::endl;
	std::cout << "min = " << highestfreq << "  " << minRes << std::endl;

	DIGFREQ2FFT_IDX(lowestfreq, len, lowIdx);
	DIGFREQ2FFT_IDX(highestfreq, len, highIdx);

	MultidimArray< std::complex<double> > fftVaux;

	MultidimArray<double> filteredMap, filtered_aux;
	filtered_aux.resizeNoCopy(resVol());
	filteredMap.initZeros(resVol());

	std::cout << "freq = " << highIdx << std::endl;
	std::cout << "freq = " << lowIdx << std::endl;

	if (onlyLocalfilter== true)
	{
		std::cout << "Applying only a local filter. The map will not be sharpened" << std::endl;

		double lastRes=1e38;

		for (double idx = lowIdx; idx < highIdx; idx++)
		{
			FFT_IDX2DIGFREQ(idx, len, freq);

			if (lastRes - sampling/freq < freq_step)
				continue;

			std::cout << "freq = " << freq << std::endl;
			std::cout << "res = " << sampling/freq << std::endl;

			freqL = freq - 0.02;
			if (freqL < 0)
				freqL = 0.001;

			freqH = freq + 0.02;

			if (freqH > 0.5)
				freqH = 0.5;

			fftVaux.initZeros(fftV);

			// Filter the input volume
			long n=0;
			double ideltal=PI/(freq-freqH);
			double ideltah=PI/(freqL-freq);


			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fftVaux)
			{
				double iun=DIRECT_MULTIDIM_ELEM(iu,n);
				double un=1.0/iun;
				if (un <= freqH && freq <= un)
				{
					//double H=0.5*(1+cos((un-w1)*ideltal));
					DIRECT_MULTIDIM_ELEM(fftVaux, n) = DIRECT_MULTIDIM_ELEM(fftV, n);
					DIRECT_MULTIDIM_ELEM(fftVaux, n) *= 0.5*(1+cos((un-freq)*ideltal));//H;
				}
				if (freqL<=un && un<=freq)
				{
					//double H=0.5*(1+cos((un-w1)*ideltal));
					DIRECT_MULTIDIM_ELEM(fftVaux, n) = DIRECT_MULTIDIM_ELEM(fftV, n);
					DIRECT_MULTIDIM_ELEM(fftVaux, n) *= 0.5*(1+cos((un-freq)*ideltah));//H;
				}
				if (freqL<=un && un<=freqH)
				{
					//double H=0.5*(1+cos((un-w1)*ideltal));
					DIRECT_MULTIDIM_ELEM(fftVaux, n) = DIRECT_MULTIDIM_ELEM(fftV, n);
				}
			}

			transformer_inv.inverseFourierTransform(fftVaux, filtered_aux);

			double stdres = sampling/sqrt(sigma);

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(filtered_aux)
			{
				double weight, sumweight, digitalfreq;
				//TODO: make the converstion in the produdesideinfo
				digitalfreq = sampling/DIRECT_MULTIDIM_ELEM(resVol, n);

				weight = exp(-(digitalfreq - freq)*(digitalfreq - freq)/(stdres*stdres));
				DIRECT_MULTIDIM_ELEM(filteredMap, n) += weight*DIRECT_MULTIDIM_ELEM(filtered_aux, n);
				sumweight += weight;
			}
			lastRes = sampling/freq;
		}

		Image<double> saveMap;
		saveMap = filteredMap;
		saveMap.write(fnOut);
	}
	else
	{
//		std::cout << "Applying only a local sharpening" << std::endl;
//
//        MultidimArray<double> auxVol;
//        MultidimArray<double> operatedfiltered, Vk, filteredVol;
//        double lastnorm=0, lastporc=1;
//        double freq;
//        double step = 0.2;
//        int idx, bool1=1, bool2=1;
//        int lastidx = -1;
//
//		size_t Niter = 5;
//
//		for (size_t i = 1; i<=Niter; ++i)
//        {
//			//std::cout << "----------------Iteration " << i << "----------------" << std::endl;
//			auxVol = filteredVol;
//			FourierTransformer transformer;
//			transformer.FourierTransform(auxVol, fftV);
//
//			localfilteredMap(fftV, operatedfiltered, minRes, maxRes, step);
//
//			filteredVol = Vorig;
//
//			filteredVol -= operatedfiltered;
//
//			//calculate norm for Vorig
//			if (i==1)
//			{
//				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vorig)
//				{
//					   normOrig +=(DIRECT_MULTIDIM_ELEM(Vorig,n)*DIRECT_MULTIDIM_ELEM(Vorig,n));
//				}
//				normOrig = sqrt(normOrig);
//				//std::cout << "norma del original  " << normOrig << std::endl;
//			}
//
//
//			//calculate norm for operatedfiltered
//			double norm=0;
//			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(operatedfiltered)
//			{
//					norm +=(DIRECT_MULTIDIM_ELEM(operatedfiltered,n)*DIRECT_MULTIDIM_ELEM(operatedfiltered,n));
//			}
//			norm=sqrt(norm);
//
//
//			double porc=lastnorm*100/norm;
//			//std::cout << "norm " << norm << " percetage " << porc << std::endl;
//
//			double subst=porc-lastporc;
//
//			if ((subst<1)&&(bool1==1)&&(i>2))
//			{
//				bool1=2;
//				//std::cout << "-----iteration completed-----" << std::endl;
//			}
//
//			lastnorm=norm;
//			lastporc=porc;
//
//			if (i==1 && lambda==1)
//			{
//				lambda=(normOrig/norm)/12;
//				std::cout << "  lambda  " << lambda << std::endl;
//			}
//
//			////Second operator
//			transformer.FourierTransform(filteredVol, fftV);
//			localfiltering(fftV, filteredVol, minRes, maxRes, step);
//
//			if (i == 1)
//					Vk = Vorig;
//			else
//					Vk = sharpenedMap;
//
//			//sharpenedMap=Vk+lambda*(filteredVol);
//			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sharpenedMap)
//			{
//				DIRECT_MULTIDIM_ELEM(sharpenedMap,n)=DIRECT_MULTIDIM_ELEM(Vk,n)+
//									 lambda*DIRECT_MULTIDIM_ELEM(filteredVol,n);
//									 //-0.01*DIRECT_MULTIDIM_ELEM(Vk,n)*SGN(DIRECT_MULTIDIM_ELEM(Vk,n));
//				if (DIRECT_MULTIDIM_ELEM(sharpenedMap,n)<-4*desvOutside_Vorig)
//					DIRECT_MULTIDIM_ELEM(sharpenedMap,n)=-4*desvOutside_Vorig;
//			}
//
////        		double desv_sharp=0;
////                computeAvgStdev_within_binary_mask(resVol, sharpenedMap, desv_sharp);
////                std::cout << "desv_sharp = " << desv_sharp << std::endl;
//
//			filteredVol = sharpenedMap;
//
//			if (bool1==2)
//			{
//				Image<double> filteredvolume;
//				filteredvolume() = sharpenedMap;
//				filteredvolume.write(fnOut);
//				break;
//			}
//        }
	}


}
