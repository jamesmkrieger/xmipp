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

#include "angular_discrete_assign2.h"
#include <data/mask.h>
#include <data/symmetries.h>
#include <reconstruction/directions.h>

// Empty constructor =======================================================
ProgAngularDiscreteAssign2::ProgAngularDiscreteAssign2()
{
    produces_a_metadata = true;
    each_image_produces_an_output = true;
    comparator = NULL;
    ctfImage = NULL;
    rank = 0;
}

ProgAngularDiscreteAssign2::~ProgAngularDiscreteAssign2()
{
	delete comparator;
	delete ctfImage;
}

// Read arguments ==========================================================
void ProgAngularDiscreteAssign2::readParams()
{
	XmippMetadataProgram::readParams();
    fnVol = getParam("--ref");
    fnSym = getParam("--sym");
    shiftStep = getDoubleParam("--shift_step");
    angStep = getDoubleParam("--ang_step");
    minTilt = getDoubleParam("--min_tilt");
    maxTilt = getDoubleParam("--max_tilt");
    maxShift = getDoubleParam("--max_shift");
    maxResol = getDoubleParam("--max_resolution");
    Ts = getDoubleParam("--sampling");
    Rmax = getIntParam("--Rmax");
    pad = getIntParam("--padding");
    ignoreCTF = checkParam("--ignoreCTF");
    phaseFlipped = checkParam("--phaseFlipped");
}

// Show ====================================================================
void ProgAngularDiscreteAssign2::show()
{
    if (!verbose)
        return;
	XmippMetadataProgram::show();
    std::cout
    << "Reference volume:    " << fnVol              << std::endl
    << "Symmetry:            " << fnSym              << std::endl
    << "Shift step:          " << shiftStep          << std::endl
    << "Angular step:        " << angStep            << std::endl
    << "Min. Tilt:           " << minTilt            << std::endl
    << "Max. Tilt:           " << maxTilt            << std::endl
    << "Max. Shift:          " << maxShift           << std::endl
    << "Max. Resolution:     " << maxResol           << std::endl
    << "Sampling:            " << Ts                 << std::endl
    << "Max. Radius:         " << Rmax               << std::endl
    << "Padding factor:      " << pad                << std::endl
    << "Ignore CTF:          " << ignoreCTF          << std::endl
    << "Phase flipped:       " << phaseFlipped       << std::endl
    ;
}

// usage ===================================================================
void ProgAngularDiscreteAssign2::defineParams()
{
    addUsageLine("Make a discrete angular assignment");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --ref <volume>              : Reference volume");
    addParamsLine("  [--shift_step <s=1>]         : Shift step");
    addParamsLine("  [--ang_step <a=5>]           : Angular step");
    addParamsLine("  [--min_tilt <t=0>]           : Minimum tilt");
    addParamsLine("  [--max_tilt <t=180>]         : Maximum tilt");
    addParamsLine("  [--max_shift <s=10>]         : Maximum shift allowed in pixels");
    addParamsLine("  [--max_resolution <f=4>]     : Maximum resolution (A)");
    addParamsLine("  [--sampling <Ts=1>]          : Sampling rate (A/pixel)");
    addParamsLine("  [--Rmax <R=-1>]              : Maximum radius (px). -1=Half of volume size");
    addParamsLine("  [--padding <p=2>]            : Padding factor");
    addParamsLine("  [--ignoreCTF]                : Ignore CTF");
    addParamsLine("  [--phaseFlipped]             : Input images have been phase flipped");
    addParamsLine("  [--sym <sym_file=c1>]        : It is used for computing the asymmetric unit");
    addExampleLine("A typical use is:",false);
    addExampleLine("xmipp_angular_discrete_assign2 -i images.xmd --ref reference.vol -o assigned_angles.xmd");
}

// Produce side information ================================================
void ProgAngularDiscreteAssign2::preProcess()
{
    // Read the reference volume
	Image<double> V;
	if (rank==0)
	{
		V.read(fnVol);
		V().setXmippOrigin();
	    Xdim=XSIZE(V());
	}
	else
	{
		size_t ydim, zdim, ndim;
		getImageSize(fnVol, Xdim, ydim, zdim, ndim);
	}

    // Construct mask
    if (Rmax<0)
    	Rmax=Xdim/2;
    Mask mask;
    mask.type = BINARY_CIRCULAR_MASK;
    mask.mode = INNER_MASK;
    mask.R1 = Rmax;
    mask.generate_mask(Xdim,Xdim);
    mask2D=mask.get_binary_mask();
    iMask2Dsum=1.0/mask2D.sum();

    // Get the list of rot, tilt, psi and shift
    SymList SL;
    if (fnSym != "")
        SL.readSymmetryFile(fnSym);
    std::vector<double> rotListAux, tiltListAux;
    make_even_distribution(rotListAux, tiltListAux, angStep, SL, false);
    size_t Nrl=rotListAux.size();
    for (size_t n=0; n<Nrl; ++n)
    {
    	if (tiltListAux[n]>=minTilt && tiltListAux[n]<=maxTilt)
    	{
    		rotList.push_back(rotListAux[n]);
    		tiltList.push_back(tiltListAux[n]);
    	}
    }

    for (double shift=-maxShift; shift<=maxShift; shift+=shiftStep)
    	shiftList.push_back(shift);

    for (double psi=0; psi<360; psi+=angStep)
    	psiList.push_back(psi);

    Ncandidates = Nrl * psiList.size() * shiftList.size() * shiftList.size();
    currentL2.initZeros(Ncandidates);
    fullL2.initZeros(Ncandidates);

    // Construct comparator
    if (rank==0)
    	comparator = new FourierComparator(V(),pad,Ts/maxResol,NEAREST); // BSPLINE3, LINEAR
    else
    	comparator = new FourierComparator(pad,Ts/maxResol,NEAREST);
}

void ProgAngularDiscreteAssign2::updateCTFImage(double defocusU, double defocusV, double angle)
{
	ctf.K=1; // get pure CTF with no envelope
	currentDefocusU=ctf.DeltafU=defocusU;
	currentDefocusV=ctf.DeltafV=defocusV;
	currentAngle=ctf.azimuthal_angle=angle;
	ctf.produceSideInfo();
	if (ctfImage==NULL)
	{
		ctfImage = new MultidimArray<double>();
		ctfImage->resizeNoCopy(comparator->volumeSize, comparator->volumeSize);
		STARTINGY(*ctfImage)=STARTINGX(*ctfImage)=0;
	}
	ctf.generateCTF(comparator->volumeSize, comparator->volumeSize, *ctfImage,Ts);
	if (phaseFlipped)
		FOR_ALL_ELEMENTS_IN_ARRAY2D(*ctfImage)
			A2D_ELEM(*ctfImage,i,j)=fabs(A2D_ELEM(*ctfImage,i,j));
}

// Predict =================================================================
void ProgAngularDiscreteAssign2::evaluateImage(const MultidimArray< std::complex<double> > &FIexp,
		                                       int idx0, int idxF)
{
	size_t NRotTilt=rotList.size();
	size_t NPsi=psiList.size();
	size_t Nshift=shiftList.size();

	size_t Nshift2=Nshift*Nshift;
	size_t Nshift2Psi=Nshift2*NPsi;

	double bestRot, bestTilt, bestPsi, bestSx, bestSy;
	double bestL2=1e38;
	size_t bestIdx=0;
//	std::cout << "New range " << idx0 << " " << idxF << std::endl;

	for (size_t nrt=0; nrt<NRotTilt; ++nrt)
	{
		size_t irt=nrt*Nshift2Psi;
		double rot=rotList[nrt];
		double tilt=tiltList[nrt];
		for (size_t npsi=0; npsi<NPsi; ++npsi)
		{
			double psi=psiList[npsi];
			comparator->setEuler(rot,tilt,psi);
			size_t irtp = irt+npsi*Nshift2;
			for (size_t nsy=0; nsy<Nshift; ++nsy)
			{
				double sy=shiftList[nsy];
				size_t irtpy = irtp+nsy*Nshift;
				for (size_t nsx=0; nsx<Nshift; ++nsx)
				{
					size_t irtpxy = irtpy+nsx;
					if (DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)>=0 && !DIRECT_MULTIDIM_ELEM(fullL2,irtpxy))
					{
						double sx=shiftList[nsx];
						double newL2=comparator->compare(FIexp,sx,sy,idx0,idxF,ctfImage);
//						double newL2=comparator->compare(FIexp,sx,sy,0,100,ctfImage);
						DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)+=newL2;
						if (DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)<bestL2)
						{
							bestRot=rot;
							bestTilt=tilt;
							bestPsi=psi;
							bestSx=sx;
							bestSy=sy;
							bestIdx=irtpxy;
							bestL2=DIRECT_MULTIDIM_ELEM(currentL2,irtpxy);
//							std::cout << rot << " " << tilt << " " << psi << " " << sy << " " << sx << " -> " << bestL2 << std::endl;
						}
					}
				}
			}
		}
	}

	// Evaluate full
	comparator->setEuler(bestRot,bestTilt,bestPsi);
	DIRECT_MULTIDIM_ELEM(currentL2,bestIdx)=comparator->compare(FIexp,bestSx,bestSy,3,comparator->maxIdx,ctfImage);
	DIRECT_MULTIDIM_ELEM(fullL2,bestIdx)=1;
//	std::cout << "Full evaluation of best " << DIRECT_MULTIDIM_ELEM(currentL2,bestIdx) << std::endl;
	bestL2 = DIRECT_MULTIDIM_ELEM(currentL2,bestIdx);

	// Remove larger items from the list of candidates
	for (size_t nrt=0; nrt<NRotTilt; ++nrt)
	{
		size_t irt=nrt*Nshift2Psi;
		double rot=rotList[nrt];
		double tilt=tiltList[nrt];
		for (size_t npsi=0; npsi<NPsi; ++npsi)
		{
			double psi=psiList[npsi];
			size_t irtp = irt+npsi*Nshift2;
			for (size_t nsy=0; nsy<Nshift; ++nsy)
			{
				double sy=shiftList[nsy];
				size_t irtpy = irtp+nsy*Nshift;
				for (size_t nsx=0; nsx<Nshift; ++nsx)
				{
					size_t irtpxy = irtpy+nsx;
					if (DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)>bestL2)
						DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)=-1;
				}
			}
		}
	}
}

//#define DEBUG
void ProgAngularDiscreteAssign2::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
    rowOut=rowIn;

    // Read input image and initial parameters
	if ((rowIn.containsLabel(MDL_CTF_DEFOCUSU) || rowIn.containsLabel(MDL_CTF_MODEL)) && !ignoreCTF)
	{
		hasCTF=true;
		ctf.readFromMdRow(rowIn);
		ctf.produceSideInfo();
	}
	else
		hasCTF=false;

	if (verbose>=2)
		std::cout << "Processing " << fnImg << std::endl;
	I.read(fnImg);

	// Apply mask
	MultidimArray<double> &mI=I();
	mI.setXmippOrigin();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mask2D)
	if (!DIRECT_MULTIDIM_ELEM(mask2D,n))
		DIRECT_MULTIDIM_ELEM(mI,n)=0.0;

	// Compute Fourier transform
	transformer2D.FourierTransform(mI, FI, false);

	try
	{
		currentL2.initZeros();
		fullL2.initZeros();
		int idxStep=4;
		for (int idx=3; idx<Xdim; idx+=idxStep)
			evaluateImage(FI, idx,idx+idxStep);
	}
	catch (XmippError XE)
	{
		std::cerr << XE << std::endl;
		std::cerr << "Warning: Cannot refine " << fnImg << std::endl;
		rowOut.setValue(MDL_ENABLED,-1);
	}
/*
	rowOut.setValue(MDL_IMAGE_ORIGINAL, fnImg);
    rowOut.setValue(MDL_IMAGE, fnImgOut);
    rowOut.setValue(MDL_ANGLE_ROT,  old_rot+p(7));
    rowOut.setValue(MDL_ANGLE_TILT, old_tilt+p(8));
    rowOut.setValue(MDL_ANGLE_PSI,  old_psi+p(9));
    rowOut.setValue(MDL_SHIFT_X,    0.);
    rowOut.setValue(MDL_SHIFT_Y,    0.);
    rowOut.setValue(MDL_FLIP,       false);
    rowOut.setValue(MDL_COST,       cost);
    if (optimizeGrayValues)
    {
		rowOut.setValue(MDL_CONTINUOUS_GRAY_A,p(0));
		rowOut.setValue(MDL_CONTINUOUS_GRAY_B,p(1));
    }
    rowOut.setValue(MDL_CONTINUOUS_SCALE_X,p(4));
    rowOut.setValue(MDL_CONTINUOUS_SCALE_Y,p(5));
    rowOut.setValue(MDL_CONTINUOUS_SCALE_ANGLE,p(6));
    rowOut.setValue(MDL_CONTINUOUS_X,p(2)+old_shiftX);
    rowOut.setValue(MDL_CONTINUOUS_Y,p(3)+old_shiftY);
    rowOut.setValue(MDL_CONTINUOUS_FLIP,old_flip);
    if (hasCTF)
    {
    	rowOut.setValue(MDL_CTF_DEFOCUSU,old_defocusU+p(10));
    	rowOut.setValue(MDL_CTF_DEFOCUSV,old_defocusV+p(11));
    	rowOut.setValue(MDL_CTF_DEFOCUS_ANGLE,old_defocusAngle+p(12));
    	if (old_defocusU+p(10)<0 || old_defocusU+p(11)<0)
    		rowOut.setValue(MDL_ENABLED,-1);
    }
    //Saving correlation and imed values in the metadata
    rowOut.setValue(MDL_CORRELATION_IDX, corrIdx);
    rowOut.setValue(MDL_CORRELATION_MASK, corrMask);
    rowOut.setValue(MDL_CORRELATION_WEIGHT, corrWeight);
    rowOut.setValue(MDL_IMED, imedDist);
*/

#ifdef DEBUG
    std::cout << "p=" << p << std::endl;
    MetaData MDaux;
    MDaux.addRow(rowOut);
    MDaux.write("PPPmd.xmd");
    Image<double> save;
    save()=P();
    save.write("PPPprojection.xmp");
    save()=I();
    save.write("PPPexperimental.xmp");
    //save()=C;
    //save.write("PPPC.xmp");
    Ip.write("PPPexperimentalp.xmp");
    Ifiltered.write("PPPexperimentalFiltered.xmp");
    Ifilteredp.write("PPPexperimentalFilteredp.xmp");
    E.write("PPPresidual.xmp");
    std::cout << A << std::endl;
    std::cout << fnImgOut << " rewritten\n";
    std::cout << "Press any key" << std::endl;
    char c; std::cin >> c;
#endif
}
#undef DEBUG

void ProgAngularDiscreteAssign2::postProcess()
{
	MetaData &ptrMdOut=*getOutputMd();
	ptrMdOut.removeDisabled();
	ptrMdOut.write(fn_out.replaceExtension("xmd"));
}
