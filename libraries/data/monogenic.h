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
#ifndef _MONOGENIC_HH
#define _MONOGENIC_HH

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <core/metadata.h>
#include <core/xmipp_image.h>
#include <data/sampling.h>
#include <core/xmipp_fft.h>
#include <core/xmipp_fftw.h>
#include <math.h>
#include <limits>
#include <complex>
#include <data/fourier_filter.h>
#include <data/filters.h>
//@{
/** Routines for working with monogenic signals
*/

class Monogenic
{
public:
	//This function determines a multidimArray with the frequency values in Fourier space;
	MultidimArray< std::complex<double> > fourierFreqs_3D(const MultidimArray< std::complex<double> > &myfftV,
			const MultidimArray<double> &inputVol);

	//This function determines a multidimArray with the frequency values in Fourier space;
	MultidimArray< std::complex<double> > fourierFreqs_2D(const MultidimArray< std::complex<double> > &myfftImg,
			const MultidimArray<double> &inputImg);

	//It computes the monogenic amplitude of an input volume;
	void monogenicAmplitude_3D(const MultidimArray<double> &inputVol, MultidimArray<double> &amplitude, int numberOfThreads);

	//Fast method: It computes the monogenic amplitude of an input volume;
	void monogenicAmplitude_3D(FourierTransformer &transformer, MultidimArray< std::complex<double> > &iu,
			Matrix1D<double> &freq_fourier_z, Matrix1D<double> &freq_fourier_y, Matrix1D<double> &freq_fourier_x,
			MultidimArray< std::complex<double> > &fftVRiesz, MultidimArray< std::complex<double> > &fftVRiesz_aux,
			MultidimArray<double> &VRiesz, const MultidimArray<double> &inputVol,
			MultidimArray< std::complex<double> > &myfftV, MultidimArray<double> &amplitude);

	//It computes the monogenic amplitude of an input image;
	void monogenicAmplitude_2D(const MultidimArray<double> &inputVol, MultidimArray<double> &amplitude, int numberOfThreads);

	//Fast method: It computes the monogenic amplitude of an input image;
	void monogenicAmplitude_2D(FourierTransformer &transformer, MultidimArray< std::complex<double> > &iu,
			Matrix1D<double> &freq_fourier_y, Matrix1D<double> &freq_fourier_x,
			MultidimArray< std::complex<double> > &fftVRiesz, MultidimArray< std::complex<double> > &fftVRiesz_aux,
			MultidimArray<double> &VRiesz, const MultidimArray<double> &inputVol,
			MultidimArray< std::complex<double> > &myfftV, MultidimArray<double> &amplitude);
};

//@}
#endif
