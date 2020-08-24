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
#include <core/xmipp_image_generic.h>
#include <data/mask.h>
#include <data/symmetries.h>
#include <reconstruction/directions.h>

// Empty constructor =======================================================
ProgAngularDiscreteAssign2::ProgAngularDiscreteAssign2()
{
    produces_a_metadata = true;
    each_image_produces_an_output = true;
    produces_an_output = true;
    output_is_stack = true;
    comparator = NULL;
    ctfImage = NULL;
    rank = 0;
}

ProgAngularDiscreteAssign2::~ProgAngularDiscreteAssign2()
{
	delete comparator;
	delete ctfImage;
	for (size_t i=0; i<shiftPhase.size(); i++)
		delete shiftPhase[i];
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
    saveReprojection = checkParam("--saveReprojection");
    saveResiduals = checkParam("--saveResidual");
    onlyEvaluate = checkParam("--onlyEvaluate");
    adjustProfile = checkParam("--adjustProfile");
    if (checkParam("--profile"))
    	fnProfile = getParam("--profile");
    if (adjustProfile)
    {
    	fnProfile = getParam("--adjustProfile",0);
    	Nadjust = getIntParam("--adjustProfile",1);

        produces_a_metadata = false;
        each_image_produces_an_output = false;
        produces_an_output = false;
	}
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
	<< "Save reprojection:   " << saveReprojection   << std::endl
	<< "Save residuals:      " << saveResiduals      << std::endl
	<< "Only evaluate:       " << onlyEvaluate       << std::endl
	<< "Adjust profile:      " << adjustProfile      << std::endl
    ;
    if (adjustProfile)
		std::cout << "Profile:             " << fnProfile << std::endl
		          << "Nadjust:             " << Nadjust   << std::endl;
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
    addParamsLine("  [--padding <p=1>]            : Padding factor");
    addParamsLine("  [--ignoreCTF]                : Ignore CTF");
    addParamsLine("  [--phaseFlipped]             : Input images have been phase flipped");
    addParamsLine("  [--saveReprojection]         : Save reprojection");
    addParamsLine("  [--saveResidual]             : Save residual");
    addParamsLine("  [--sym <sym_file=c1>]        : It is used for computing the asymmetric unit");
    addParamsLine("  [--onlyEvaluate]             : Only evaluate, assumes that there are angles and shifts");
    addParamsLine("  [--profile <profile>]        : Adjusted profile to apply to the volume");
    addParamsLine("  [--adjustProfile <profile=\"\"> <N=1000>] : Adjust volume Fourier amplitude with the first N images");
    addExampleLine("A typical use is:",false);
    addExampleLine("xmipp_angular_discrete_assign2 -i images.xmd --ref reference.vol -o assigned_angles.xmd");
}

// Produce side information ================================================
void ProgAngularDiscreteAssign2::preProcess()
{
    // Read the reference volume
	Image<double> V;
	V.read(fnVol);
	V().setXmippOrigin();
	Xdim=XSIZE(V());

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
    maskbg.initZeros(mask2D);
    stdBgSum=stdBgN=0;
    stdFgSum=stdFgN=0;

    if (!onlyEvaluate)
    {
		// Get the list of rot, tilt, psi and shift
		SymList SL;
		if (fnSym != "")
			SL.readSymmetryFile(fnSym);
		std::vector<double> rotListAux, tiltListAux;
		make_even_distribution(rotListAux, tiltListAux, angStep, SL, true);
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

		if (fnProfile!="")
		{
			Matrix1D<double> vProfile;
			vProfile.readNoSize(fnProfile);
			profile=vProfile;
		}
    }

    // Construct comparator
    comparator = new FourierComparator(V(),pad,Ts/maxResol,NEAREST); // BSPLINE3, LINEAR
    comparator->KV=profile;

    // Precompute phase planes
    if (!onlyEvaluate)
    {
		for (size_t nsy=0; nsy<shiftList.size(); ++nsy)
		{
			double sy=shiftList[nsy];
			for (size_t nsx=0; nsx<shiftList.size(); ++nsx)
			{
				double sx=shiftList[nsx];
				MultidimArray< std::complex<double> > *phasePlane=new MultidimArray< std::complex<double> >();
				comparator->preparePhasePlane(sx, sy, *phasePlane);
				shiftPhase.push_back(phasePlane);
			}
		}
    }
    else
    {
    	MultidimArray< std::complex<double> > *phasePlane=new MultidimArray< std::complex<double> >();
    	shiftPhase.push_back(phasePlane);
    }

	// Create output files
	fnMask = fn_out.removeAllExtensions()+"_mask.stk";
	fnWeight = fn_out.removeAllExtensions()+"_fourierMask.stk";
	if (saveReprojection)
		fnReprojection = fn_out.removeAllExtensions()+"_reprojections.stk";
	if (saveResiduals)
		fnResidual = fn_out.removeAllExtensions()+"_residuals.stk";
	if (rank==0 && !onlyEvaluate)
	{
		createEmptyFile(fnWeight, XSIZE(comparator->weightFourier), YSIZE(comparator->weightFourier), 1, mdInSize, true, WRITE_OVERWRITE);
		createEmptyFile(fnMask, Xdim, Xdim, 1, mdInSize, true, WRITE_OVERWRITE);
		if (saveReprojection)
			createEmptyFile(fnReprojection, Xdim, Xdim, 1, mdInSize, true, WRITE_OVERWRITE);
		if (saveResiduals)
			createEmptyFile(fnResidual, Xdim, Xdim, 1, mdInSize, true, WRITE_OVERWRITE);
	}
    filter.FilterShape = REALGAUSSIAN;
    filter.FilterBand = LOWPASS;
    filter.w1 = 2;
}

void ProgAngularDiscreteAssign2::updateCTFImage()
{
	ctf.produceSideInfo();
	if (ctfImage==NULL)
	{
		ctfImage = new MultidimArray<double>();
		ctfImage->resizeNoCopy(comparator->volumeSize, comparator->volumeSize);
		STARTINGY(*ctfImage)=STARTINGX(*ctfImage)=0;
	}
	ctf.generateCTF(comparator->volumeSize, comparator->volumeSize, *ctfImage, Ts);
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

	double bestL2=1e38;
	size_t bestIdx=0;
	bestIdxPhase=0;
	//std::cout << "New range " << idx0 << " " << idxF << std::endl;
	comparator->prepareForIdx(idx0,idxF);

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
			size_t idxPhasePlane=0;
			for (size_t nsy=0; nsy<Nshift; ++nsy)
			{
				double sy=shiftList[nsy];
				size_t irtpy = irtp+nsy*Nshift;
				for (size_t nsx=0; nsx<Nshift; ++nsx, ++idxPhasePlane)
				{
					size_t irtpxy = irtpy+nsx;
					if (DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)>=0 && !DIRECT_MULTIDIM_ELEM(fullL2,irtpxy))
					{
						double sx=shiftList[nsx];
						double newL2=comparator->compare(FIexp,shiftPhase[idxPhasePlane],ctfImage);
						DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)+=newL2;
						if (DIRECT_MULTIDIM_ELEM(currentL2,irtpxy)<bestL2)
						{
							bestRot=rot;
							bestTilt=tilt;
							bestPsi=psi;
							bestSx=sx;
							bestSy=sy;
							bestIdx=irtpxy;
							bestIdxPhase=idxPhasePlane;
							bestL2=DIRECT_MULTIDIM_ELEM(currentL2,irtpxy);
							//std::cout << rot << " " << tilt << " " << psi << " " << sy << " " << sx << " -> " << bestL2 << std::endl;
						}
					}
				}
			}
		}
	}

	// Evaluate full
	comparator->setEuler(bestRot,bestTilt,bestPsi);
	comparator->prepareForIdx(3,comparator->maxIdx);
	DIRECT_MULTIDIM_ELEM(currentL2,bestIdx)=comparator->compare(FIexp,shiftPhase[bestIdxPhase],ctfImage);
	DIRECT_MULTIDIM_ELEM(fullL2,bestIdx)=1;
	//std::cout << "Full evaluation of best " << DIRECT_MULTIDIM_ELEM(currentL2,bestIdx) << std::endl;
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

void ProgAngularDiscreteAssign2::evaluateResiduals(const MultidimArray< std::complex<double> > &FIexp)
{
	comparator->prepareForIdx(0,XSIZE(FIexp));
	comparator->setEuler(bestRot,bestTilt,bestPsi);
	comparator->compare(FI,shiftPhase[bestIdxPhase],ctfImage,true);

	wFI()=comparator->weightFourier;
	P()=comparator->projection();
	const MultidimArray<double> &mP=P();
	double th=mP.computeMax()/20;
	maskbg.initZeros();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mP)
	if (DIRECT_MULTIDIM_ELEM(mask2D,n)==0 || DIRECT_MULTIDIM_ELEM(mP,n)<th)
		DIRECT_MULTIDIM_ELEM(maskbg,n)=1;

	FI=comparator->shiftedExpFourier;
	transformer2D.inverseFourierTransform(); // Substitute I by its shifted version
	double avg, stddev;
	I().computeAvgStdev_within_binary_mask(maskbg, avg, stddev);
	stdBgSum+=stddev;
	stdBgN+=1;

	stddev=stdBgSum/stdBgN;
	maskDampen().resizeNoCopy(maskbg);
	maskfg.initZeros(maskbg);
	maskDampen().initConstant(1.0);
	const MultidimArray<double> &mD=maskDampen();
	const MultidimArray<double> &mI=I();
	Qbg=0;
	double Nbg=0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mP)
	{
		if (DIRECT_MULTIDIM_ELEM(maskbg,n)==1 || DIRECT_MULTIDIM_ELEM(mask2D,n)==0)
		{
			double z=DIRECT_MULTIDIM_ELEM(mI,n)/stddev;
			DIRECT_MULTIDIM_ELEM(mD,n)=exp(-0.5*z*z);
			Qbg+=DIRECT_MULTIDIM_ELEM(mD,n);
			Nbg+=1;
		}
		else
			DIRECT_MULTIDIM_ELEM(maskfg,n)=1;
	}
	QCCreal=correlationIndex(mP,mI,&mask2D);
	QFourier=wFI().computeAvg();
	Qbg/=Nbg;

	typeCast(maskfg,maskfgD);
    filter.applyMaskSpace(maskfgD);

	P().rangeAdjustLS(I(),&maskfg);
	double Nfg=0;
	Qfg=0;
	IP()=mI;
	const MultidimArray<double> &mIP=IP();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mP)
	{
		double diff=DIRECT_MULTIDIM_ELEM(mI,n)-DIRECT_MULTIDIM_ELEM(mP,n);
		DIRECT_MULTIDIM_ELEM(mIP,n)-=DIRECT_MULTIDIM_ELEM(maskfgD,n)*DIRECT_MULTIDIM_ELEM(mP,n);

		if (DIRECT_MULTIDIM_ELEM(maskfg,n)==1)
		{
			double absdiff=std::abs(diff);
			if (stdFgN<1e6)
			{
				stdFgSum+=absdiff;
				stdFgN+=1;
			}

			double z=absdiff/(3*stdFgSum/stdFgN);
			Qfg+=exp(-0.5*z*z);
			Nfg+=1;
		}

	}
	Qfg/=Nfg;
}

void ProgAngularDiscreteAssign2::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut)
{
	if (adjustProfile)
	{
		if (Nadjust==0)
			return;
		else
			Nadjust--;
	}

    rowOut=rowIn;

    // Read input image and initial parameters
	if ((rowIn.containsLabel(MDL_CTF_DEFOCUSU) || rowIn.containsLabel(MDL_CTF_MODEL)) && !ignoreCTF)
	{
		hasCTF=true;
		ctf.readFromMdRow(rowIn);
		ctf.produceSideInfo();
		updateCTFImage();
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
		if (!onlyEvaluate)
		{
			currentL2.initZeros();
			fullL2.initZeros();
			int idxStep=4;
			for (int idx=3; idx<Xdim; idx+=idxStep)
				evaluateImage(FI, idx,idx+idxStep);

			rowOut.setValue(MDL_ANGLE_ROT,  bestRot);
			rowOut.setValue(MDL_ANGLE_TILT, bestTilt);
			rowOut.setValue(MDL_ANGLE_PSI,  bestPsi);
			rowOut.setValue(MDL_SHIFT_X,    bestSx);
			rowOut.setValue(MDL_SHIFT_Y,    bestSy);
		}
		else
		{
			rowIn.getValue(MDL_ANGLE_ROT,  bestRot);
			rowIn.getValue(MDL_ANGLE_TILT, bestTilt);
			rowIn.getValue(MDL_ANGLE_PSI,  bestPsi);
			rowIn.getValue(MDL_SHIFT_X,    bestSx);
			rowIn.getValue(MDL_SHIFT_Y,    bestSy);

			comparator->preparePhasePlane(bestSx, bestSy, *(shiftPhase[0]));
			bestIdxPhase=0;
		}
		evaluateResiduals(FI);
		if (!onlyEvaluate)
	    	rowOut.setValue(MDL_MAXCC, QCCreal);
    	rowOut.setValue(MDL_QCCREAL, QCCreal);
    	rowOut.setValue(MDL_QBACKGROUND, Qbg);
    	rowOut.setValue(MDL_QFOREGROUND, Qfg);
    	rowOut.setValue(MDL_QFOURIER, QFourier);

		// Write Fourier mask
		size_t idx=fnImgOut.getPrefixNumber();
		FileName fn;
		fn.compose(idx,fnWeight);
		wFI.write(fn,0,true,WRITE_OVERWRITE);
		rowOut.setValue(MDL_IMAGE_MASK_FOURIER, fn);

		fn.compose(idx,fnMask);
		maskDampen.write(fn,0,true,WRITE_OVERWRITE);
		rowOut.setValue(MDL_IMAGE_MASK, fn);

	    // Write reprojection
	    if (saveReprojection)
		{
			FileName fn;
			fn.compose(idx,fnReprojection);
			P.write(fn,0,true,WRITE_OVERWRITE);
		    rowOut.setValue(MDL_IMAGE_REF, fn);
		}
	    if (saveResiduals)
	    {
			FileName fn;
	    	fn.compose(idx,fnResidual);
	    	IP.write(fn,0,true,WRITE_OVERWRITE);
		    rowOut.setValue(MDL_IMAGE_RESIDUAL, fn);
	    }
	}
	catch (XmippError XE)
	{
		std::cerr << XE << std::endl;
		std::cerr << "Warning: Cannot refine " << fnImg << std::endl;
		rowOut.setValue(MDL_ENABLED,-1);
	}
}

void ProgAngularDiscreteAssign2::finishProcessing()
{
	if (adjustProfile)
	{
		profile.initZeros(XSIZE(comparator->IabsSum));
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(comparator->IabsSum)
			if (DIRECT_MULTIDIM_ELEM(comparator->VabsSum,n)>0)
				DIRECT_MULTIDIM_ELEM(profile,n)=DIRECT_MULTIDIM_ELEM(comparator->IabsSum,n)/DIRECT_MULTIDIM_ELEM(comparator->VabsSum,n);
	}
}

void ProgAngularDiscreteAssign2::postProcess()
{
	if (!adjustProfile)
	{
		MetaData &ptrMdOut=*getOutputMd();
		ptrMdOut.removeDisabled();
		ptrMdOut.write(fn_out.replaceExtension("xmd"));
	}
	else
		profile.write(fnProfile);
}
