/***************************************************************************
 *
 * Authors:    Jose Luis Vilas, 					  jlvilas@cnb.csic.es
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

#ifndef _PROG_RESBFACTOR
#define _PROG_RESBFACTOR

#include <iostream>
#include <core/xmipp_program.h>
#include <core/xmipp_image.h>
#include <core/metadata.h>
#include <core/xmipp_fft.h>
#include <core/xmipp_fftw.h>
#include <math.h>
#include <limits>
#include <complex>
#include <data/fourier_filter.h>
#include <data/filters.h>
#include <string>


/**@defgroup Resolution B-Factor
   @ingroup ReconsLibrary */
//@{

class ProgResBFactor : public XmippProgram
{
public:
	 /** Filenames */
	FileName fnOut, fn_pdb, fn_locres;

	/** sampling rate*/
	double sampling;

	/** Number of atoms in the pdb or alpha-carbons*/
	int numberOfAtoms;

	bool medianTrue;

	std::vector<double> residuesToChimera;
	double fscResolution;

	struct pdbdata
	{
		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> z;
		std::vector<double> b;
		std::vector<int> residue;
		std::vector<double> atomCovRad;
	};
	struct pdbdata at_pos;

public:

    void defineParams();
    void readParams();
    void produceSideInfo();
    void analyzePDB();

    template<typename T>
    std::vector<size_t> sort_indexes(const std::vector<T> &v);

    void sweepByResidue(MultidimArray<int> &mask, std::vector<double> &residuesToChimera);

    template<typename T>
    void maskFromPDBData(struct pdbdata &coord, MultidimArray<T> &mask);

    void generateOutputPDB(std::vector<double> &residuesToChimera);

    void smoothingResidueOutput();

    void run();
};
//@}
#endif
