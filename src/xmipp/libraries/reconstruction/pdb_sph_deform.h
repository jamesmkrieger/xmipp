/***************************************************************************
 *
 * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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
#ifndef _PROG_PDB_SPH_DEFORM
#define _PROG_PDB_SPH_DEFORM

#include <data/pdb.h>
#include <core/xmipp_program.h>


class ProgPdbSphDeform: public XmippProgram
{
public:
    /** PDB file */
    FileName fn_pdb;

    /** Deformation coefficients list */
    FileName fn_sph;

    /** Output fileroot */
    FileName fn_out;

    //Vector containing the deformation coefficients
	Matrix1D<double> clnm;

public:
    /** Params definitions */
    void defineParams();

    /** Read from a command line. */
    void readParams();

    /** Show parameters. */
    void show();

    /** Run. */
    void run();

    /** Number of coefficients */
    int numberCoefficients();
};
//@}
#endif