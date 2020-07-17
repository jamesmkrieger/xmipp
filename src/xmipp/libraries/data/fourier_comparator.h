/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
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
 * 02111-1307  USAcd ..
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef _CORE_FOURIER_FOURIER_COMPARATOR_H
#define _CORE_FOURIER_FOURIER_COMPARATOR_H

#include <core/xmipp_image.h>
#include <core/xmipp_fftw.h>

/**@defgroup FourierComparator Fourier comparator
   @ingroup ReconsLibrary */
//@{

/// @defgroup Projections Comparator in Fourier space
/// @ingroup DataLibrary
//@{
/** Program class to create projections in Fourier space */
class FourierComparator
{
public:
    /// Padding factor
    double paddingFactor;
    /// Maximum Frequency for pixels
    double maxFrequency;
    /// The order of B-Spline for interpolation
    double BSplineDeg;
public:
    // Volume to project
    MultidimArray<double> *volume;

    // Real and imaginary B-spline coefficients for Fourier of the volume
    MultidimArray< double > VfourierRealCoefs, VfourierImagCoefs;

    // Phase shift image
    MultidimArray<double> phaseShiftImgB, phaseShiftImgA;

    // Original volume size
    int volumeSize;

    // Volume padded size
    int volumePaddedSize;

    // Euler matrix
    Matrix2D<float> E;

    // Fourier index
    MultidimArray<int> wIdx;
    MultidimArray<float> wx, wy;
    MultidimArray<int> wIdxXcount, wIdxYcount;
    int idx0, idxF;

    // MaxIndex depending on maxFrequency
    int maxIdx;

    // Transformer
    FourierTransformer transformer2D;

    // Projection in Fourier space
    MultidimArray< std::complex<double> > projectionFourier;

    // Projection
    Image<double> projection;
public:
    /* Empty constructor */
    FourierComparator(double paddFactor, double maxFreq, int degree);

    /*
     * The constructor of the class
     */
    FourierComparator(MultidimArray<double> &V, double paddFactor, double maxFreq, int degree);

	/** Set Euler matrix */
    void setEuler(double rot, double tilt, double psi);

    /** Prepare phase plane for a shift */
    void preparePhasePlane(double shiftx, double shifty, MultidimArray< std::complex<double> > &phase);

    /** Prepare phase plane for a shift */
    void prepareForIdx(int idx0, int idxF);

    /**
     * Compare
     */
    double compare(const MultidimArray< std::complex<double> > &Iexp,
                   const MultidimArray< std::complex<double> > *phaseShift=NULL,
				   const MultidimArray<double> *ctf=NULL);

    /** Update volume */
    void updateVolume(MultidimArray<double> &V);
public:
    /// Prepare the Spline coefficients and projection space
    void produceSideInfo();

    /// Prepare projection space
    void produceSideInfoProjection();
};
//@}

#endif
