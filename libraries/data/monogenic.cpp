/***************************************************************************
 *
 * Authors:     Jose Luis Vilas (jlvilas@cnb.csic.es)
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
#include "monogenic.h"


//This function determines a multidimArray with the frequency values in Fourier space;
MultidimArray< std::complex<double> > Monogenic::fourierFreqs_3D(const MultidimArray< std::complex<double> > &myfftV,
		const MultidimArray<double> &inputVol)
{
	MultidimArray< std::complex<double> > iu;

	iu.initZeros(myfftV);

	double uz, uy, ux, uz2, u2, uz2y2;
	long n=0;

	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
		uz2=uz*uz;

		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
			uz2y2=uz2+uy*uy;

			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				FFT_IDX2DIGFREQ(j,XSIZE(inputVol), ux);
				u2=uz2y2+ux*ux;
				if ((k != 0) || (i != 0) || (j != 0))
					DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
				else
					DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
				++n;
			}
		}
	}

	return iu;
}

MultidimArray< std::complex<double> > Monogenic::fourierFreqs_2D(const MultidimArray< std::complex<double> > &myfftImg,
		const MultidimArray<double> &inputImg)
{
	MultidimArray< std::complex<double> > iu;

	iu.initZeros(myfftImg);

	double uz, uy, ux, uy2, u2;
	long n=0;

	for(size_t i=0; i<YSIZE(myfftImg); ++i)
	{
		FFT_IDX2DIGFREQ(i,YSIZE(inputImg),uy);
		uy2 = uy*uy;

		for(size_t j=0; j<XSIZE(myfftImg); ++j)
		{
			FFT_IDX2DIGFREQ(j,XSIZE(inputImg), ux);
			u2=uy2+ux*ux;
			if ((i != 0) || (j != 0))
				DIRECT_MULTIDIM_ELEM(iu,n) = 1.0/sqrt(u2);
			else
				DIRECT_MULTIDIM_ELEM(iu,n) = 1e38;
			++n;
		}
	}

	return iu;
}

void Monogenic::monogenicAmplitude_3D(const MultidimArray<double> &inputVol, MultidimArray<double> &amplitude, int numberOfThreads)
{
	MultidimArray< std::complex<double> > iu, fftVRiesz, fftVRiesz_aux, myfftV;
	MultidimArray<double> VRiesz;
	VRiesz= inputVol;



	FourierTransformer transformer;
	transformer.setThreadsNumber(numberOfThreads);

	transformer.FourierTransform(VRiesz, myfftV);

	double u;
	Matrix1D<double> freq_fourier_z, freq_fourier_y, freq_fourier_x;

	freq_fourier_z.initZeros(ZSIZE(myfftV));
	freq_fourier_x.initZeros(XSIZE(myfftV));
	freq_fourier_y.initZeros(YSIZE(myfftV));

	VEC_ELEM(freq_fourier_z,0) = 1e-38;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol), u);
		VEC_ELEM(freq_fourier_z,k) = u;
	}
	VEC_ELEM(freq_fourier_y,0) = 1e-38;
	for(size_t k=0; k<YSIZE(myfftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,YSIZE(inputVol), u);
		VEC_ELEM(freq_fourier_y,k) = u;
	}
	VEC_ELEM(freq_fourier_x,0) = 1e-38;
	for(size_t k=0; k<XSIZE(myfftV); ++k)
	{
		FFT_IDX2DIGFREQ(k,XSIZE(inputVol), u);
		VEC_ELEM(freq_fourier_x,k) = u;
	}

	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);

	iu = fourierFreqs_3D(myfftV, inputVol);

	monogenicAmplitude_3D(transformer, iu, freq_fourier_z, freq_fourier_y, freq_fourier_x,
							fftVRiesz, fftVRiesz_aux, VRiesz, inputVol, myfftV, amplitude);
}

void Monogenic::monogenicAmplitude_3D(FourierTransformer &transformer, MultidimArray< std::complex<double> > &iu,
		Matrix1D<double> &freq_fourier_z, Matrix1D<double> &freq_fourier_y, Matrix1D<double> &freq_fourier_x,
		MultidimArray< std::complex<double> > &fftVRiesz, MultidimArray< std::complex<double> > &fftVRiesz_aux,
		MultidimArray<double> &VRiesz, const MultidimArray<double> &inputVol,
		MultidimArray< std::complex<double> > &myfftV, MultidimArray<double> &amplitude)
{
	VRiesz.resizeNoCopy(inputVol);

	std::complex<double> J(0,1);

	fftVRiesz.initZeros(myfftV);
	fftVRiesz_aux.initZeros(myfftV);

	//Original volume in real space and preparing Riesz components
	long n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(iu,n);
				++n;
			}
		}
	}

	transformer.inverseFourierTransform(fftVRiesz, amplitude);

//	#ifdef DEBUG_DIR
//		Image<double> filteredvolume;
//		filteredvolume = VRiesz;
//		filteredvolume.write(formatString("Volumen_filtrado_%i_%i.vol", dir+1,count));
//	#endif

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n) *= DIRECT_MULTIDIM_ELEM(amplitude,n);


	// Calculate first component of Riesz vector
	double ux;
	n=0;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				ux = VEC_ELEM(freq_fourier_x,j);
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}

	transformer.inverseFourierTransform(fftVRiesz, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	// Calculate second and third component of Riesz vector
	n=0;
	double uy, uz;
	for(size_t k=0; k<ZSIZE(myfftV); ++k)
	{
		uz = VEC_ELEM(freq_fourier_z,k);
		for(size_t i=0; i<YSIZE(myfftV); ++i)
		{
			uy = VEC_ELEM(freq_fourier_y,i);
			for(size_t j=0; j<XSIZE(myfftV); ++j)
			{
				DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = uz*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
				++n;
			}
		}
	}

	transformer.inverseFourierTransform(fftVRiesz, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n) += DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	transformer.inverseFourierTransform(fftVRiesz_aux, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n) += DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
}


//void Monogenic::monogenicAmplitudeHPF_3D(const MultidimArray<double> &inputVol, MultidimArray<double> &amplitude, int numberOfThreads=1)
//{
//	MultidimArray< std::complex<double> > iu, fftVRiesz, fftVRiesz_aux, myfftV;
//	MultidimArray<double> VRiesz;
//	VRiesz.resizeNoCopy(inputVol);
//
//	FourierTransformer transformer;
//	transformer.setThreadsNumber(numberOfThreads);
//
//	transformer.FourierTransform(inputVol, myfftV);
//
//	double u;
//	Matrix1D<double> freq_fourier_z, freq_fourier_y, freq_fourier_x;
//
//	freq_fourier_z.initZeros(ZSIZE(myfftV));
//	freq_fourier_x.initZeros(XSIZE(myfftV));
//	freq_fourier_y.initZeros(YSIZE(myfftV));
//
//	VEC_ELEM(freq_fourier_z,0) = 1e-38;
//	for(size_t k=0; k<ZSIZE(myfftV); ++k)
//	{
//		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol), u);
//		VEC_ELEM(freq_fourier_z,k) = u;
//	}
//	VEC_ELEM(freq_fourier_y,0) = 1e-38;
//	for(size_t k=0; k<YSIZE(myfftV); ++k)
//	{
//		FFT_IDX2DIGFREQ(k,YSIZE(inputVol), u);
//		VEC_ELEM(freq_fourier_y,k) = u;
//	}
//	VEC_ELEM(freq_fourier_x,0) = 1e-38;
//	for(size_t k=0; k<XSIZE(myfftV); ++k)
//	{
//		FFT_IDX2DIGFREQ(k,XSIZE(inputVol), u);
//		VEC_ELEM(freq_fourier_x,k) = u;
//	}
//
//	fftVRiesz.initZeros(myfftV);
//	fftVRiesz_aux.initZeros(myfftV);
//
//	iu = fourierFreqs_3D(myfftV, inputVol);
//
//	monogenicAmplitude_3D(transformer, iu, freq_fourier_z, freq_fourier_y, freq_fourier_x,
//							fftVRiesz, fftVRiesz_aux, VRiesz, inputVol, myfftV, amplitude);
//}



void Monogenic::monogenicAmplitude_2D(const MultidimArray<double> &inputImg, MultidimArray<double> &amplitude, int numberOfThreads)
{
	MultidimArray< std::complex<double> > iu, fftVRiesz, fftVRiesz_aux, myfftImg;
	MultidimArray<double> VRiesz;
	VRiesz = inputImg;

	FourierTransformer transformer;
	transformer.setThreadsNumber(numberOfThreads);

	transformer.FourierTransform(VRiesz, myfftImg);

	double u;
	Matrix1D<double> freq_fourier_y, freq_fourier_x;

	freq_fourier_x.initZeros(XSIZE(myfftImg));
	freq_fourier_y.initZeros(YSIZE(myfftImg));

	VEC_ELEM(freq_fourier_y,0) = 1e-38;
	for(size_t k=0; k<YSIZE(myfftImg); ++k)
	{
		FFT_IDX2DIGFREQ(k,YSIZE(inputImg), u);
		VEC_ELEM(freq_fourier_y,k) = u;
	}
	VEC_ELEM(freq_fourier_x,0) = 1e-38;
	for(size_t k=0; k<XSIZE(myfftImg); ++k)
	{
		FFT_IDX2DIGFREQ(k,XSIZE(inputImg), u);
		VEC_ELEM(freq_fourier_x,k) = u;
	}

	fftVRiesz.initZeros(myfftImg);
	fftVRiesz_aux.initZeros(myfftImg);

	iu = fourierFreqs_2D(myfftImg, inputImg);

	monogenicAmplitude_2D(transformer, iu, freq_fourier_y, freq_fourier_x,
							fftVRiesz, fftVRiesz_aux, VRiesz, inputImg, myfftImg, amplitude);
}

void Monogenic::monogenicAmplitude_2D(FourierTransformer &transformer, MultidimArray< std::complex<double> > &iu,
		Matrix1D<double> &freq_fourier_y, Matrix1D<double> &freq_fourier_x,
		MultidimArray< std::complex<double> > &fftVRiesz, MultidimArray< std::complex<double> > &fftVRiesz_aux,
		MultidimArray<double> &VRiesz, const MultidimArray<double> &inputImg,
		MultidimArray< std::complex<double> > &myfftImg, MultidimArray<double> &amplitude)
{
	VRiesz.resizeNoCopy(inputImg);

	std::complex<double> J(0,1);

	fftVRiesz.initZeros(myfftImg);
	fftVRiesz_aux.initZeros(myfftImg);

	//Original volume in real space and preparing Riesz components
	long n=0;
	for(size_t i=0; i<YSIZE(myfftImg); ++i)
	{
		for(size_t j=0; j<XSIZE(myfftImg); ++j)
		{
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = DIRECT_MULTIDIM_ELEM(myfftImg, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = -J;
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(fftVRiesz, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) *= DIRECT_MULTIDIM_ELEM(iu,n);
			++n;
		}
	}

	transformer.inverseFourierTransform(fftVRiesz, amplitude);

//	#ifdef DEBUG_DIR
//		Image<double> filteredvolume;
//		filteredvolume = VRiesz;
//		filteredvolume.write(formatString("Volumen_filtrado_%i_%i.vol", dir+1,count));
//	#endif

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n) *= DIRECT_MULTIDIM_ELEM(amplitude,n);


	// Calculate first and second component of Riesz vector
	double ux, uy;
	n=0;
	for(size_t i=0; i<YSIZE(myfftImg); ++i)
	{
		uy = VEC_ELEM(freq_fourier_y,i);
		for(size_t j=0; j<XSIZE(myfftImg); ++j)
		{
			ux = VEC_ELEM(freq_fourier_x,j);
			DIRECT_MULTIDIM_ELEM(fftVRiesz, n) = ux*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
			DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n) = uy*DIRECT_MULTIDIM_ELEM(fftVRiesz_aux, n);
			++n;
		}
	}

	transformer.inverseFourierTransform(fftVRiesz, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n)+=DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);

	transformer.inverseFourierTransform(fftVRiesz_aux, VRiesz);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(amplitude)
		DIRECT_MULTIDIM_ELEM(amplitude,n) += DIRECT_MULTIDIM_ELEM(VRiesz,n)*DIRECT_MULTIDIM_ELEM(VRiesz,n);
		DIRECT_MULTIDIM_ELEM(amplitude,n)=sqrt(DIRECT_MULTIDIM_ELEM(amplitude,n));
}

//void Monogenic::bandPassFilterFunction(FourierTransformer &transformer, MultidimArray< std::complex<double> > &iu,
//		const MultidimArray< std::complex<double> > &myfftV,
//                double w, double wL, MultidimArray<double> &filteredVol, int count)
//{
//		MultidimArray< std::complex<double> > &fftVfilter;
//        fftVfilter.initZeros(myfftV);
//        size_t xdim, ydim, zdim, ndim;
//
//        double delta = wL-w;
//        double w_inf = w-delta;
//        // Filter the input volume and add it to amplitude
//        long n=0;
//        double ideltal=PI/(delta);
//        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(myfftV)
//        {
//                double un=DIRECT_MULTIDIM_ELEM(iu,n);
//                if (un>=w && un<=wL)
//                {
//                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
//                } else if (un<=w && un>=w_inf)
//                {
//                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) = DIRECT_MULTIDIM_ELEM(myfftV, n);
//                        DIRECT_MULTIDIM_ELEM(fftVfilter, n) *= 0.5*(1+cos((un-w)*ideltal));//H;
//                }
//        }
//
//        filteredVol.resizeNoCopy(Vorig);
//
//        transformer.inverseFourierTransform(fftVfilter, filteredVol);
//
////        #ifdef DEBUG_FILTER
////        Image<double> filteredvolume;
////        filteredvolume() = filteredVol;
////        filteredvolume.write(formatString("Volumen_filtrado_%i.vol", count));
////        #endif
//}
