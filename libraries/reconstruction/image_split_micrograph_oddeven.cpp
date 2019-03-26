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

#include "image_split_micrograph_oddeven.h"
#include "classify_extract_features.h"
#include <data/filters.h>
#include <data/normalize.h>
#include <core/histogram.h>
#include <fstream>

// Read arguments ==========================================================
void ProgSplitMic::readParams()
{
    fnOdd = getParam("--odd");
    fnEven = getParam("--even");
    fnOut = getParam("-o");

}

// Show ====================================================================
void ProgSplitMic::showInfo()
{
    if (verbose==0)
        return;
    std::cerr
    << "Odd Image:           " << fnOdd         << std::endl
    << "Even Image:          " << fnEven        << std::endl
    << "Output Image:        " << fnOut       << std::endl;
}

// Usage ===================================================================
void ProgSplitMic::defineParams()
{
    addUsageLine("Given two aligned micrograph from frames it cast a reliability image");
    addParamsLine("  --odd <micrograph>                       : odd micrograph");
    addParamsLine("  --even <micrograph>                      : even micrograph");
    addParamsLine("  -o <micrograph>                          : Output - reliability image");
}

//void ProgWeightCdf::monogenicAmplitude2D()
//{
//
//
//
//}



void ProgSplitMic::run()
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

	noiseImage = oddImage();
	noiseImage -= evenImage();
	avgImage = oddImage();
	avgImage += evenImage();
	normalize_OldXmipp(noiseImage);
	normalize_OldXmipp(avgImage);
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


}
