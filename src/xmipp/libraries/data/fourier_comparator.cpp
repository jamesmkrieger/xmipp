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
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "fourier_comparator.h"
#include <core/xmipp_fft.h>
#include <core/geometry.h>
#include <core/transformations.h>
#include <core/bilib/kernel.h>

FourierComparator::FourierComparator(double paddFactor, double maxFreq, int degree)
{
    paddingFactor = paddFactor;
    maxFrequency = maxFreq;
    BSplineDeg = degree;
    volume = NULL;
}

FourierComparator::FourierComparator(MultidimArray<double> &V, double paddFactor, double maxFreq, int degree)
{
    paddingFactor = paddFactor;
    maxFrequency = maxFreq;
    BSplineDeg = degree;
    updateVolume(V);
}

void FourierComparator::updateVolume(MultidimArray<double> &V)
{
    volume = &V;
    volumeSize=XSIZE(*volume);
    produceSideInfo();
}

void FourierComparator::setEuler(double rot, double tilt, double psi)
{
	Matrix2D<double> Ed;
    Euler_angles2matrix(rot,tilt,psi,Ed);
    typeCast(Ed,E);
}

void FourierComparator::preparePhasePlane(double shiftx, double shifty, MultidimArray< std::complex<double> > &phase)
{
	phase.resizeNoCopy(wIdx);
    for (size_t i=0; i<YSIZE(wIdx); ++i)
    {
        double yyshifti = -2 * PI * shifty / volumeSize * i;
        for (size_t j=0; j<XSIZE(wIdx); ++j)
        {
            double xxshift = -2 * PI * shiftx / volumeSize;

            double a,b;
            sincos(xxshift*j+yyshifti,&b,&a);
			DIRECT_A2D_ELEM(phase,i,j) = std::complex<double>(a,b);
        }
    }
}

/** Prepare phase plane for a shift */
void FourierComparator::prepareForIdx(int iidx0, int iidxF)
{
	idx0=iidx0;
	idxF=iidxF;

	wIdxYcount.initZeros(YSIZE(wIdx));
	wIdxXcount.initZeros(XSIZE(wIdx));

	int idxFF=std::min(idxF,maxIdx);
    for (size_t i=0; i<YSIZE(wIdx); ++i)
    {
        for (size_t j=0; j<idxFF; ++j)
        {
        	int idxij = DIRECT_A2D_ELEM(wIdx,i,j);
        	if (idxij<idx0 || idxij>idxFF)
        		continue;
        	DIRECT_A1D_ELEM(wIdxYcount,i)++;
        	DIRECT_A1D_ELEM(wIdxXcount,j)++;
        }
    }
}

//#define DEBUG
double FourierComparator::compare(const MultidimArray< std::complex<double> > &Iexp,
								  const MultidimArray< std::complex<double> > *phaseShift,
								  const MultidimArray<double> *ctf,
								  bool keepImages)
{
	if (idx0>maxIdx)
		return 0.0;

	if (keepImages)
	{
		projectionFourier.initZeros();
		shiftedExpFourier.initZeros();
		weightFourier.initConstant(1.0);
	}

    int Xdim=(int)XSIZE(VfourierRealCoefs);
    int Ydim=(int)YSIZE(VfourierRealCoefs);
    int Zdim=(int)ZSIZE(VfourierRealCoefs);

    int idxFF=std::min(idxF,maxIdx);
    double retval = 0.0;
    for (size_t i=0; i<YSIZE(wIdx); ++i)
    {
    	if (DIRECT_A1D_ELEM(wIdxYcount,i)==0)
    		continue;

    	float freqy = DIRECT_A1D_ELEM(wy,i);

    	float freqYvol_X=MAT_ELEM(E,1,0)*freqy;
    	float freqYvol_Y=MAT_ELEM(E,1,1)*freqy;
    	float freqYvol_Z=MAT_ELEM(E,1,2)*freqy;
        for (size_t j=0; j<idxFF; ++j)
        {
        	int idxij = DIRECT_A2D_ELEM(wIdx,i,j);
        	if (idxij<idx0 || idxij>idxFF)
        		continue;

        	if (DIRECT_A1D_ELEM(wIdxXcount,j)==0)
        		continue;

            // The frequency of pairs (i,j) in 2D
        	float freqx = DIRECT_A1D_ELEM(wx,j);

            // Compute corresponding frequency in the volume
        	float freqvol_X=freqYvol_X+MAT_ELEM(E,0,0)*freqx;
        	float freqvol_Y=freqYvol_Y+MAT_ELEM(E,0,1)*freqx;
        	float freqvol_Z=freqYvol_Z+MAT_ELEM(E,0,2)*freqx;

            double c,d;

            int kVolume, iVolume, jVolume;
            if (BSplineDeg==0)
            {
                // 0 order interpolation
                // Compute corresponding index in the volume
                kVolume=(int)std::round(freqvol_Z*volumePaddedSize);
                iVolume=(int)std::round(freqvol_Y*volumePaddedSize);
                jVolume=(int)std::round(freqvol_X*volumePaddedSize);
                c = A3D_ELEM(VfourierRealCoefs,kVolume,iVolume,jVolume);
                d = A3D_ELEM(VfourierImagCoefs,kVolume,iVolume,jVolume);
            }
            else if (BSplineDeg==1)
            {
                // B-spline linear interpolation
                double kVolumed=freqvol_Z*volumePaddedSize;
                double iVolumed=freqvol_Y*volumePaddedSize;
                double jVolumed=freqvol_X*volumePaddedSize;
                c=VfourierRealCoefs.interpolatedElement3D(jVolumed,iVolumed,kVolumed);
                d=VfourierImagCoefs.interpolatedElement3D(jVolumed,iVolumed,kVolumed);
                kVolume=(int)std::round(kVolumed);
                iVolume=(int)std::round(iVolumed);
                jVolume=(int)std::round(jVolumed);
            }
            else
            {
                // B-spline cubic interpolation
                double kVolumed=freqvol_Z*volumePaddedSize;
                double iVolumed=freqvol_Y*volumePaddedSize;
                double jVolumed=freqvol_X*volumePaddedSize;

                // Commented for speed-up, the corresponding code is below
                // c=VfourierRealCoefs.interpolatedElementBSpline3D(jVolume,iVolume,kVolume);
                // d=VfourierImagCoefs.interpolatedElementBSpline3D(jVolume,iVolume,kVolume);

                // The code below is a replicate for speed reasons of interpolatedElementBSpline3D
                double z=kVolumed;
                double y=iVolumed;
                double x=jVolumed;

                // Logical to physical
                z -= STARTINGZ(VfourierRealCoefs);
                y -= STARTINGY(VfourierRealCoefs);
                x -= STARTINGX(VfourierRealCoefs);

                int l1 = (int)ceil(x - 2);
                int l2 = l1 + 3;

                int m1 = (int)ceil(y - 2);
                int m2 = m1 + 3;

                int n1 = (int)ceil(z - 2);
                int n2 = n1 + 3;

                c = d = 0.0;
                double aux;
                for (int nn = n1; nn <= n2; nn++)
                {
                    int equivalent_nn=nn;
                    if      (nn<0)
                        equivalent_nn=-nn-1;
                    else if (nn>=Zdim)
                        equivalent_nn=2*Zdim-nn-1;
                    double yxsumRe = 0.0, yxsumIm = 0.0;
                    for (int m = m1; m <= m2; m++)
                    {
                        int equivalent_m=m;
                        if      (m<0)
                            equivalent_m=-m-1;
                        else if (m>=Ydim)
                            equivalent_m=2*Ydim-m-1;
                        double xsumRe = 0.0, xsumIm = 0.0;
                        for (int l = l1; l <= l2; l++)
                        {
                            double xminusl = x - (double) l;
                            int equivalent_l=l;
                            if      (l<0)
                                equivalent_l=-l-1;
                            else if (l>=Xdim)
                                equivalent_l=2*Xdim-l-1;
                            double CoeffRe = (double) DIRECT_A3D_ELEM(VfourierRealCoefs,equivalent_nn,equivalent_m,equivalent_l);
                            double CoeffIm = (double) DIRECT_A3D_ELEM(VfourierImagCoefs,equivalent_nn,equivalent_m,equivalent_l);
                            BSPLINE03(aux,xminusl);
                            xsumRe += CoeffRe * aux;
                            xsumIm += CoeffIm * aux;
                        }

                        double yminusm = y - (double) m;
                        BSPLINE03(aux,yminusm);
						yxsumRe += xsumRe * aux;
						yxsumIm += xsumIm * aux;
                    }

                    double zminusn = z - (double) nn;
                    BSPLINE03(aux,zminusn);
					c += yxsumRe * aux;
					d += yxsumIm * aux;
                }
                kVolume=(int)std::round(kVolumed);
                iVolume=(int)std::round(iVolumed);
                jVolume=(int)std::round(jVolumed);
            }

            // Phase shift to move the origin of the image to the corner
            double a=DIRECT_A2D_ELEM(phaseShiftImgA,i,j);
            double b=DIRECT_A2D_ELEM(phaseShiftImgB,i,j);
            if (ctf!=NULL)
            {
            	double ctfij=DIRECT_A2D_ELEM(*ctf,i,j);
            	a*=ctfij;
            	b*=ctfij;
            }

            // Multiply Fourier coefficient in volume times phase shift
            double ac = a * c;
            double bd = b * d;
            double ab_cd = (a + b) * (c + d);
            double reTh = ac - bd;
            double imTh = ab_cd - ac - bd;
            if (XSIZE(KV)>0)
            {
            	double K=DIRECT_A1D_ELEM(KV,idxij);
            	reTh*=K;
            	imTh*=K;
            }

			double *ptrI_ij=(double *)&DIRECT_A2D_ELEM(Iexp,i,j);
            double reExp, imExp;
            if (phaseShift!=NULL)
            {
				c = *ptrI_ij;
				d = *(ptrI_ij+1);

				double *ptrPhase_ij=(double *)&DIRECT_A2D_ELEM(*phaseShift,i,j);
				a=*ptrPhase_ij;
				b=*(ptrPhase_ij+1);
				ac = a * c;
				bd = b * d;
				ab_cd = (a + b) * (c + d);
				reExp = ac - bd;
				imExp = ab_cd - ac - bd;
            }
            else
            {
				reExp = *ptrI_ij;
				imExp = *(ptrI_ij+1);
            }

            double reDiff = reTh-reExp;
            double imDiff = imTh-imExp;
            retval+=reDiff*reDiff+imDiff*imDiff;
            if (keepImages)
            {
            	DIRECT_A2D_ELEM(projectionFourier,i,j)=std::complex<double>(reTh,imTh);
            	DIRECT_A2D_ELEM(shiftedExpFourier,i,j)=std::complex<double>(reExp,imExp);

            	int idx=DIRECT_A2D_ELEM(wIdx,i,j);
            	DIRECT_A1D_ELEM(VabsSum,idx)+=std::abs(DIRECT_A2D_ELEM(projectionFourier,i,j));
            	DIRECT_A1D_ELEM(IabsSum,idx)+=std::abs(DIRECT_A2D_ELEM(shiftedExpFourier,i,j));

            	A3D_ELEM(VN,kVolume,iVolume,jVolume)++;
            	A3D_ELEM(VreSum,kVolume,iVolume,jVolume)+=reExp;
            	A3D_ELEM(VimSum,kVolume,iVolume,jVolume)+=imExp;
            	A3D_ELEM(VreSum2,kVolume,iVolume,jVolume)+=reExp*reExp;
            	A3D_ELEM(VimSum2,kVolume,iVolume,jVolume)+=imExp*imExp;
            	A3D_ELEM(VreimSum,kVolume,iVolume,jVolume)+=reExp*imExp;

            	if (A3D_ELEM(VN,kVolume,iVolume,jVolume)>10 && i!=0 && j!=0)
            	{
            		double iN=1.0/A3D_ELEM(VN,kVolume,iVolume,jVolume);
					double reMean = A3D_ELEM(VreSum,kVolume,iVolume,jVolume)*iN;
					double imMean = A3D_ELEM(VimSum,kVolume,iVolume,jVolume)*iN;
					double reSigma2 = A3D_ELEM(VreSum2,kVolume,iVolume,jVolume)*iN-reMean*reMean;
					double imSigma2 = A3D_ELEM(VimSum2,kVolume,iVolume,jVolume)*iN-imMean*imMean;
					double reimCov = A3D_ELEM(VreimSum,kVolume,iVolume,jVolume)*iN-reMean*imMean;

					double A=reSigma2;
					double B=reimCov;
					double D=imSigma2;
					double det=A*D-B*B;
					double X=reExp-reMean;
					double Y=imExp-imMean;

					double mahalanobisDistance = (D*X*X-2*B*X*Y+A*Y*Y)/det;
					if (mahalanobisDistance>0)
					{
						mahalanobisDistance = sqrt(mahalanobisDistance);
						if (mahalanobisDistance>2.8) // 2.8 = sqrt(2+3*sqrt(2*p))
						{
							DIRECT_A2D_ELEM(weightFourier,i,j)=exp(-mahalanobisDistance*mahalanobisDistance/(2*0.94*1.5*0.94*1.5)); // 0.94 = sqrt(2+3*sqrt(2*p))/3 for p = 2
																																	// 1.5 is a safety factor
//							std::cout << "N=" << A3D_ELEM(VN,kVolume,iVolume,jVolume) << " re=" << reExp << " im=" << imExp << " reMean=" << reMean << " imMean=" << imMean << std::endl
//									  << "   re2=" << reSigma2 << " im2=" << imSigma2 << " cov=" << reimCov << " i=" << i << " j=" << j << " dist=" << mahalanobisDistance
//									  << " weight=" << DIRECT_A2D_ELEM(weightFourier,i,j) << std::endl;
						}
					}
            	}
            }
        }
    }

    if (keepImages)
    {
	    transformer2D.inverseFourierTransform();
#ifdef DEBUG
	    // projection.write("PPPth.xmp");
		// projection.write("PPPExp.xmp");
		std::cout << "Press any key" << std::endl;
		Image<double> save;
		save()=VN; save.write("PPPVN.vol");
		save()=VreSum; save.write("PPPVreSum.vol");
		save()=VimSum; save.write("PPPVimSum.vol");
		save()=VreSum2; save.write("PPPVreSum2.vol");
		save()=VimSum2; save.write("PPPVimSum2.vol");
		save()=VreimSum; save.write("PPPVreimSum.vol");
		save()=weightFourier; save.write("PPPMahalanobis.xmp");
		std::cout << "Press any key" << std::endl;
		char ch; std::cin >> ch;
#endif
	}
    return retval;
}
#undef DEBUG

void FourierComparator::produceSideInfo()
{
    // Zero padding
    MultidimArray<double> Vpadded;
    int paddedDim=(int)(paddingFactor*volumeSize);
    volume->window(Vpadded,FIRST_XMIPP_INDEX(paddedDim),FIRST_XMIPP_INDEX(paddedDim),FIRST_XMIPP_INDEX(paddedDim),
                   LAST_XMIPP_INDEX(paddedDim),LAST_XMIPP_INDEX(paddedDim),LAST_XMIPP_INDEX(paddedDim));
    volume->clear();

    // Make Fourier transform, shift the volume origin to the volume center and center it
    MultidimArray< std::complex<double> > Vfourier;
    FourierTransformer transformer3D;
    transformer3D.completeFourierTransform(Vpadded,Vfourier);
    ShiftFFT(Vfourier, FIRST_XMIPP_INDEX(XSIZE(Vpadded)), FIRST_XMIPP_INDEX(YSIZE(Vpadded)), FIRST_XMIPP_INDEX(ZSIZE(Vpadded)));
    CenterFFT(Vfourier,true);
    Vfourier.setXmippOrigin();

    // Compensate for the Fourier normalization factor
    double K=(double)(XSIZE(Vpadded)*XSIZE(Vpadded)*XSIZE(Vpadded))/(double)(volumeSize*volumeSize);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Vfourier)
    DIRECT_MULTIDIM_ELEM(Vfourier,n)*=K;
    Vpadded.clear();

    // Compute Bspline coefficients
    if (BSplineDeg==3)
    {
        MultidimArray< double > VfourierRealAux, VfourierImagAux;
        Complex2RealImag(Vfourier, VfourierRealAux, VfourierImagAux);
        Vfourier.clear();
        produceSplineCoefficients(BSPLINE3,VfourierRealCoefs,VfourierRealAux);

        // Release memory as soon as you can
        VfourierRealAux.clear();

        // Remove all those coefficients we are sure we will not use during the projections
        volumePaddedSize=XSIZE(VfourierRealCoefs);
        int idxMax=maxFrequency*XSIZE(VfourierRealCoefs)+10; // +10 is a safety guard
        idxMax=std::min(FINISHINGX(VfourierRealCoefs),idxMax);
        int idxMin=std::max(-idxMax,STARTINGX(VfourierRealCoefs));
        VfourierRealCoefs.selfWindow(idxMin,idxMin,idxMin,idxMax,idxMax,idxMax);

        produceSplineCoefficients(BSPLINE3,VfourierImagCoefs,VfourierImagAux);
        VfourierImagAux.clear();
        VfourierImagCoefs.selfWindow(idxMin,idxMin,idxMin,idxMax,idxMax,idxMax);
    }
    else {
        Complex2RealImag(Vfourier, VfourierRealCoefs, VfourierImagCoefs);
        volumePaddedSize=XSIZE(VfourierRealCoefs);
    }

    VreSum.initZeros(VfourierRealCoefs);
    VimSum.initZeros(VfourierRealCoefs);
    VreSum2.initZeros(VfourierRealCoefs);
    VimSum2.initZeros(VfourierRealCoefs);
    VreimSum.initZeros(VfourierRealCoefs);
    VN.initZeros(VfourierRealCoefs);

    produceSideInfoProjection();
}

void FourierComparator::produceSideInfoProjection()
{
    // Auxiliary FFT transformer
    projection().initZeros(volumeSize,volumeSize);
    projection().setXmippOrigin();
    transformer2D.FourierTransform(projection(),projectionFourier,false);
    shiftedExpFourier.resizeNoCopy(projectionFourier);
    weightFourier.resizeNoCopy(projectionFourier);

    // Calculate phase shift terms
    phaseShiftImgA.initZeros(projectionFourier);
    phaseShiftImgB.initZeros(projectionFourier);
    wIdx.initZeros(projectionFourier);
    double shift=-FIRST_XMIPP_INDEX(volumeSize);
    double xxshift = -2 * PI * shift / volumeSize;
    wy.initZeros(YSIZE(projectionFourier));
    wx.initZeros(XSIZE(projectionFourier));
    for (size_t i=0; i<YSIZE(projectionFourier); ++i)
    {
        double phasey=(double)(i) * xxshift;
        double wyi;
    	FFT_IDX2DIGFREQ(i,volumeSize,wyi);
    	DIRECT_A1D_ELEM(wy,i)=(float)wyi;
        for (size_t j=0; j<XSIZE(projectionFourier); ++j)
        {
            // Phase shift to move the origin of the image to the corner
            double dotp = (double)(j) * xxshift + phasey;
            sincos(dotp,&DIRECT_A2D_ELEM(phaseShiftImgB,i,j),&DIRECT_A2D_ELEM(phaseShiftImgA,i,j));

            double wxj;
        	FFT_IDX2DIGFREQ(j,volumeSize,wxj);
        	DIRECT_A1D_ELEM(wx,j)=(float)wxj;
        	A2D_ELEM(wIdx,i,j)=(int)std::round(std::sqrt(wxj*wxj+wyi*wyi)*volumeSize);
        }
    }

    maxIdx = (int)std::ceil(maxFrequency*volumeSize);
    VabsSum.initZeros(maxIdx+1);
    IabsSum.initZeros(maxIdx+1);
}
