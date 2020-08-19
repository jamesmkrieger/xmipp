/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#ifndef _PROG_ANGULAR_PREDICT_DISCRETE2
#define _PROG_ANGULAR_PREDICT_DISCRETE2

#include <core/xmipp_metadata_program.h>
#include <data/ctf.h>
#include <data/fourier_comparator.h>
#include <data/fourier_filter.h>

/**@defgroup AngularPredictDiscretes2 angular_discrete_assign2 (Discrete angular assignment2)
   @ingroup ReconsLibrary */
//@{

/** Predict Continuous Parameters. */
class ProgAngularDiscreteAssign2: public XmippMetadataProgram
{
public:
    /** Filename of the reference volume */
    FileName fnVol;
    /** Filename of the symmetry file */
    FileName fnSym;
    /** Shift step */
    double shiftStep;
    /** Angular step */
    double angStep;
    /** Minimum tilt */
    double minTilt;
    /** Maximum tilt */
    double maxTilt;
    /** Maximum shift allowed */
    double maxShift;
    /** Maximum frequency (A) */
    double maxResol;
    /** Sampling rate */
    double Ts;
    /** Maximum radius */
    int Rmax;
    /** Padding factor */
    int pad;
    // Ignore CTF
    bool ignoreCTF;
    // Phase Flipped
    bool phaseFlipped;
    // Save reprojection
    bool saveReprojection;
    // Save differences
    bool saveResiduals;
    // Only evaluate
    bool onlyEvaluate;
    // Adjust profile
    bool adjustProfile;
    // Number of images to adjust profile
    int Nadjust;
    // Output profile
    FileName fnProfile;
public:
    // Best alignment
	double bestRot, bestTilt, bestPsi, bestSx, bestSy;
	// Best phase plane
	size_t bestIdxPhase;
	// Auxiliary variable for composing output names
	FileName fnWeight, fnMask, fnReprojection, fnResidual;

	// Rank (used for MPI version)
    int rank;

    // 2D mask in real space
    MultidimArray<int> mask2D;
    // Inverse of the sum of Mask2D
    double iMask2Dsum;
    // Fourier projector
    FourierComparator *comparator;
    // Volume size
    size_t Xdim;
    // Input image
	Image<double> I;
    // Projection image
	Image<double> P;
    // IP, residual I-P
	Image<double> IP;
    // Weight in Fourier
	Image<double> wFI;
	// Mask
	MultidimArray<int> maskbg, maskfg;
	// Mask smooth borders
	MultidimArray<double> maskfgD;
	// Filter
	FourierFilter filter;
	// Dampening mask
	Image<double> maskDampen;
	// Has CTF
	bool hasCTF;
	// CTF
	CTFDescription ctf;
	// Current defoci
	double currentDefocusU, currentDefocusV, currentAngle;
	// CTF image
	MultidimArray<double> *ctfImage;
	// Projection directions
    std::vector<double> rotList, tiltList, psiList;
	// Shift list
    std::vector<double> shiftList;
    // Number of candidates
    size_t Ncandidates;
	// Current distance
    MultidimArray<double> currentL2;
    // Distance is full
    MultidimArray<unsigned char> fullL2;

    // Cross correlation
    double QCCreal, Qbg, Qfg, QFourier;

    // Std sum
    double stdBgSum, stdBgN;
    double stdFgSum, stdFgN;

    // Profile of the volume
	MultidimArray<double> profile;

    // Transformer
    FourierTransformer transformer2D;

    // Fourier transform of I
    MultidimArray< std::complex<double> > FI;

    // Shift phase planes
    std::vector< MultidimArray<std::complex<double> >* > shiftPhase;

	// Amplitude correction profile
    MultidimArray<double> amplitudeProfile;
public:
    /// Empty constructor
    ProgAngularDiscreteAssign2();

    /// Destructor
    ~ProgAngularDiscreteAssign2();

    /// Read argument from command line
    void readParams();

    /// Show
    void show();

    /// Define parameters
    void defineParams();

    /** Produce side info. All ranks.
        An exception is thrown if any of the files is not found*/
    void preProcess();

	/** Evaluate image between indexes [idx0, idxF) */
    void evaluateImage(const MultidimArray< std::complex<double> > &FIexp, int idx0, int idxF);

    /** Evaluate residuals */
    void evaluateResiduals(const MultidimArray< std::complex<double> > &FIexp);

    /** Predict angles and shift.
        At the input the pose parameters must have an initial guess of the
        parameters. At the output they have the estimated pose.*/
    void processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut);

    /** Update CTF image */
    void updateCTFImage(double defocusU, double defocusV, double angle);

    /** Finish processing */
    void finishProcessing();

    /** Post process */
    void postProcess();
};
//@}
#endif
