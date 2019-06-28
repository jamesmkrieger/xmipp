/***************************************************************************
 *
 * Authors:    Jose Luis Vilas,                     jlvilas@cnb.csic.es
 * 			   Carlos Oscar Sorzano					coss@cnb.csic.es
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

#include "volume_directional_sharpening.h"
#include "resolution_directional.h"
//#define DEBUG
//#define DEBUG_MASK

void ProgDirSharpening::readParams()
{
        fnVol = getParam("--vol");
        fnMask = getParam("--mask");
        sampling = getDoubleParam("--sampling");
        res_step = getDoubleParam("--resStep");
        significance = getDoubleParam("--significance");
        R = getDoubleParam("--volumeRadius");
        lambda = getDoubleParam("-l");
        K= getDoubleParam("-k");
        Niter = getIntParam("-i");
        Nthread = getIntParam("-n");
        fnOut = getParam("-o");
        fnMD = getParam("--md");
}

void ProgDirSharpening::defineParams()
{
        addUsageLine("This function performs local sharpening");
        addParamsLine("  --vol <vol_file=\"\">   : Input volume");
        addParamsLine("  --mask <vol_file=\"\">  : Binary mask");
        addParamsLine("  --sampling <s=1>: sampling");
        addParamsLine("  [--volumeRadius <s=100>]                : This parameter determines the radius of a sphere where the volume is");
        addParamsLine("  [--significance <s=0.95>]               : The level of confidence for the hypothesis test.");
        addParamsLine("  [--resStep <s=0.5>]  		             : Resolution step (precision) in A");
        addParamsLine("  -o <output=\"Sharpening.vol\">: sharpening volume");
        addParamsLine("  [--md <output=\"params.xmd\">]: sharpening params");
        addParamsLine("  [-l <lambda=1>]: regularization param");
        addParamsLine("  [-k <K=0.025>]: K param");
        addParamsLine("  [-i <Niter=50>]: iteration");
        addParamsLine("  [-n <Nthread=1>]: threads number");
}


void ProgDirSharpening::icosahedronVertex(Matrix2D<double> &vertex)
{
	std::cout << "Defining Icosahedron vertex..." << std::endl;

	//The icosahedron vertex are located in (0, +-1, +-phi), (+-1, +-phi, 0), (+-phi, 0, +-1) with phi = (1+sqrt(5))/2
	double phi =  (1+sqrt(5))/2;

	vertex.initZeros(12,3);

	MAT_ELEM(vertex, 0,0) = 0;    		MAT_ELEM(vertex, 0,1) = 1;    		MAT_ELEM(vertex, 0,2) = phi;
	MAT_ELEM(vertex, 1,0) = 0;    		MAT_ELEM(vertex, 1,1) = 1;    		MAT_ELEM(vertex, 1,2) = -phi;
	MAT_ELEM(vertex, 2,0) = 0;    		MAT_ELEM(vertex, 2,1) = -1;    		MAT_ELEM(vertex, 2,2) = phi;
	MAT_ELEM(vertex, 3,0) = 0;    		MAT_ELEM(vertex, 3,1) = -1;    		MAT_ELEM(vertex, 3,2) = -phi;
	MAT_ELEM(vertex, 4,0) = 1;    		MAT_ELEM(vertex, 4,1) = phi;    		MAT_ELEM(vertex, 4,2) = 0;
	MAT_ELEM(vertex, 5,0) = 1;    		MAT_ELEM(vertex, 5,1) = -phi;    		MAT_ELEM(vertex, 5,2) = 0;
	MAT_ELEM(vertex, 6,0) = -1;    		MAT_ELEM(vertex, 6,1) = phi;    		MAT_ELEM(vertex, 6,2) = 0;
	MAT_ELEM(vertex, 7,0) = -1;    		MAT_ELEM(vertex, 7,1) = -phi;    		MAT_ELEM(vertex, 7,2) = 0;
	MAT_ELEM(vertex, 8,0) = phi;    		MAT_ELEM(vertex, 8,1) = 0;    		MAT_ELEM(vertex, 8,2) = 1;
	MAT_ELEM(vertex, 9,0) = phi;    		MAT_ELEM(vertex, 9,1) = 0;    		MAT_ELEM(vertex, 9,2) = -1;
	MAT_ELEM(vertex, 10,0) = -phi;    		MAT_ELEM(vertex, 10,1) = 0;    		MAT_ELEM(vertex, 10,2) = 1;
	MAT_ELEM(vertex, 11,0) = -phi;    		MAT_ELEM(vertex, 11,1) = 0;    		MAT_ELEM(vertex, 11,2) = -1;

	vertex = vertex*(1/sqrt(1+phi*phi));
}

void ProgDirSharpening::icosahedronFaces(Matrix2D<int> &faces, Matrix2D<double> &vertex)
{
	std::cout << " Defining the faces of the icosahedron ..." << std::endl;
	//Each face is defined by three vertex

	//An icosahedron has 20 faces.
	faces.initZeros(20,3);

	int xdim = MAT_YSIZE(vertex); //Number of vertex
	int counter = 0;


	std::cout << "Number of vertex = " << xdim << std::endl;

	double norm_vertex = MAT_ELEM(vertex, 0,0)*MAT_ELEM(vertex, 0,0) + \
    		MAT_ELEM(vertex, 0,1)*MAT_ELEM(vertex, 0,1) +\
			MAT_ELEM(vertex, 0,2)*MAT_ELEM(vertex, 0,2);

	for (int i = 0; i<(xdim-2); ++i)
	{
	    for (int j = (i+1); j<(xdim-1); ++j)
	    {
	        for (int k = (j+1); k<(xdim); ++k)
	        {
	            if ((i==j) || (j ==k) || (i==k))
	                continue;

//	            std::cout << i << " " << j << " " << k << std::endl;
	            double dotprodutij, dotprodutjk, dotprodutik;
	            dotprodutij = (MAT_ELEM(vertex, i,0)*MAT_ELEM(vertex, j,0) + \
	            		MAT_ELEM(vertex, i,1)*MAT_ELEM(vertex, j,1) +\
						MAT_ELEM(vertex, i,2)*MAT_ELEM(vertex, j,2))/norm_vertex;

	            dotprodutjk = (MAT_ELEM(vertex, k,0)*MAT_ELEM(vertex, j,0) + \
	            	            		MAT_ELEM(vertex, k,1)*MAT_ELEM(vertex, j,1) + \
	            						MAT_ELEM(vertex, k,2)*MAT_ELEM(vertex, j,2))/norm_vertex;

	            dotprodutik = (MAT_ELEM(vertex, i,0)*MAT_ELEM(vertex, k,0) + \
	            	            		MAT_ELEM(vertex, i,1)*MAT_ELEM(vertex, k,1) + \
	            						MAT_ELEM(vertex, i,2)*MAT_ELEM(vertex, k,2))/norm_vertex;

	            // the number 65 comes because is greter than 60 that is the exact angle between two icosahedron vertex
	            if ((acos(dotprodutij)< 65*PI/180) && (acos(dotprodutjk)< 65*PI/180) && (acos(dotprodutik)< 65*PI/180) )
	            {
	            	MAT_ELEM(faces, counter, 0) = i;
	            	MAT_ELEM(faces, counter, 1) = j;
	            	MAT_ELEM(faces, counter, 2) = k;
	            	std::cout << i << " " << j << " " << k << std::endl;
	            	++counter;
	            }

	        }
	    }
	}
	//However, only the half of the sphere is used, so 10 faces must be considered

	int v1, v2, v3, v1_bis, v2_bis, v3_bis;
	double x1, x2, x3, y1, y2, y3, z1, z2, z3, x1_bis, x2_bis, x3_bis, y1_bis, y2_bis, y3_bis, z1_bis, z2_bis, z3_bis;

	std::vector<int> facestoSkip;

	for (int f1 = 0; f1<(MAT_YSIZE(faces)-1); ++f1)
	{
		v1 = MAT_ELEM(faces,f1, 0); v2 = MAT_ELEM(faces,f1, 1); v3 = MAT_ELEM(faces,f1, 2);

		for (int f2 = f1+1; f2<MAT_YSIZE(faces); ++f2)
		{
			v1_bis = MAT_ELEM(faces,f2, 0); v2_bis = MAT_ELEM(faces,f2, 1); v3_bis = MAT_ELEM(faces,f2, 2);

			if ( (v1 == v1_bis) || (v1 == v2_bis) || (v1 == v3_bis) || (v2 == v2_bis) || (v2 == v3_bis) || (v3 == v3_bis))
			{
				continue;
			}
			else
			{
				double x_tot = x1 + x2 + x3;
				double y_tot = y1 + y2 + y3;
				double z_tot = z1 + z2 + z3;

				double x_tot_bis = x1 + x2 + x3;
				double y_tot_bis = y1 + y2 + y3;
				double z_tot_bis = z1 + z2 + z3;

				x1 = MAT_ELEM(vertex,v1, 0); y1 = MAT_ELEM(vertex,v1, 1); z1 = MAT_ELEM(vertex,v1, 2);
				x2 = MAT_ELEM(vertex,v2, 0); y2 = MAT_ELEM(vertex,v2, 1); z2 = MAT_ELEM(vertex,v2, 2);
				x3 = MAT_ELEM(vertex,v3, 0); y3 = MAT_ELEM(vertex,v3, 1); z3 = MAT_ELEM(vertex,v3, 2);

				x1_bis = MAT_ELEM(vertex,v1, 0); y1_bis = MAT_ELEM(vertex,v1, 1); z1_bis = MAT_ELEM(vertex,v1, 2);
				x2_bis = MAT_ELEM(vertex,v2, 0); y2_bis = MAT_ELEM(vertex,v2, 1); z2_bis = MAT_ELEM(vertex,v2, 2);
				x3_bis = MAT_ELEM(vertex,v3, 0); y3_bis = MAT_ELEM(vertex,v3, 1); z3_bis = MAT_ELEM(vertex,v3, 2);

				if ( (abs(x_tot-x_tot_bis)<0.1 ) && (abs(y_tot-y_tot_bis)<0.1 ) && (abs(z_tot-z_tot_bis)<0.1 ))
				{
					MAT_ELEM(faces,f2, 0) = -1;
					MAT_ELEM(faces,f2, 1) = -1;
					MAT_ELEM(faces,f2, 2) = -1;
				}
			}

		}
	}



}

//void ProgDirSharpening::definingIcosahedronCone(MultidimArray< std::complex<double> > &myfftV,
//			MultidimArray< std::complex<double> > &conefilter, double rot, double tilt, Matrix2D<int> &vertex)
//	{
//	//	conefilter.initZeros(myfftV);
//		conefilter = myfftV;
//		// Filter the input volume and add it to amplitude
//
//	//	MultidimArray<double> conetest;
//	//	conetest.initZeros(myfftV);
//	//	#ifdef DEBUG_DIR
//	//	MultidimArray<double> coneVol;
//	//	coneVol.initZeros(iu);
//	//	#endif
//
//		int xdim = MAT_YSIZE(vertex);
//
//		for (size_t v = 0; v<xdim; ++v)
//		{
//
//		}
//
//		double x_dir, y_dir, z_dir;
//
//		x_dir = sin(tilt*PI/180)*cos(rot*PI/180);
//		y_dir = sin(tilt*PI/180)*sin(rot*PI/180);
//		z_dir = cos(tilt*PI/180);
//
//	//	double ang_con = 10*PI/180;
//		double ang_con = 15*PI/180;
//
//		double uz, uy, ux;
//		long n = 0;
//		for(size_t k=0; k<ZSIZE(myfftV); ++k)
//		{
//			uz = VEC_ELEM(freq_fourier,k);
//			uz *= z_dir;
//			for(size_t i=0; i<YSIZE(myfftV); ++i)
//			{
//				uy = VEC_ELEM(freq_fourier,i);
//				uy *= y_dir;
//				for(size_t j=0; j<XSIZE(myfftV); ++j)
//				{
//					double iun=DIRECT_MULTIDIM_ELEM(iu,n);
//					ux = VEC_ELEM(freq_fourier,j);
//					ux *= x_dir;
//
//					//BE CAREFULL with the order
//					//double dotproduct = (uy*x_dir + ux*y_dir + uz*z_dir)*iun;
//					iun *= (ux + uy + uz);
//					double acosine = acos(fabs(iun));
//	//				DIRECT_MULTIDIM_ELEM(conetest, n) = log(real(conj(DIRECT_MULTIDIM_ELEM(myfftV, n))*DIRECT_MULTIDIM_ELEM(myfftV, n)));
//					if (acosine>ang_con)
//					{
//						DIRECT_MULTIDIM_ELEM(conefilter, n) = 0;
//	//					DIRECT_MULTIDIM_ELEM(conetest, n) = 0;
//					}
//	/*
//					//4822.53 mean a smoothed cone angle of 20 degrees
//					double arg_exp = acosine*acosine*acosine*acosine*4822.53;
//					DIRECT_MULTIDIM_ELEM(conefilter, n) *= exp(-arg_exp);
//	//				DIRECT_MULTIDIM_ELEM(conetest, n) = exp(-arg_exp);
//	 */
//					++n;
//				}
//			}
//		}
//	//
//	//	Image<double> saveImg2;
//	//	saveImg2 = conetest;
//	//	saveImg2.write("cono.vol");
//
//
//}


void ProgDirSharpening::produceSideInfo()
{
        std::cout << "Starting..." << std::endl;

    	Image<double> V;
    	V.read(fnVol);
    	V().setXmippOrigin();

    	FourierTransformer transformer;

    	MultidimArray<double> &inputVol = V();
    	VRiesz.resizeNoCopy(inputVol);
    	maxRes = ZSIZE(inputVol);
    	minRes = 2*sampling;

    	//TODO: remove Nthr
    	int Nthr = 1;
    	transformer_inv.setThreadsNumber(Nthr);

    	Vorig = inputVol;

    	transformer.FourierTransform(inputVol, fftV);
    	iu.initZeros(fftV);

    	// Frequency volume
    	double uz, uy, ux, uz2, u2, uz2y2;
    	long n=0;

    	for(size_t k=0; k<ZSIZE(fftV); ++k)
    	{
    		FFT_IDX2DIGFREQ(k,ZSIZE(inputVol),uz);
    		uz2=uz*uz;

    		for(size_t i=0; i<YSIZE(fftV); ++i)
    		{
    			FFT_IDX2DIGFREQ(i,YSIZE(inputVol),uy);
    			uz2y2=uz2+uy*uy;

    			for(size_t j=0; j<XSIZE(fftV); ++j)
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

    	// Prepare mask
    	MultidimArray<int> &pMask=mask();
   		mask.read(fnMask);
    	mask().setXmippOrigin();

    	N_smoothing = 7;
    	NVoxelsOriginalMask = 0;
    	double radius = 0;
    	FOR_ALL_ELEMENTS_IN_ARRAY3D(pMask)
    	{
    		if (A3D_ELEM(pMask, k, i, j) == 1)
    		{
    			if ((k*k + i*i + j*j)>radius)
    				radius = k*k + i*i + j*j;
    		}
    //		std::cout << "i j k " << i << " " << j << " " << k << std::endl;

    		if (A3D_ELEM(pMask, k, i, j) == 1)
    			++NVoxelsOriginalMask;
    		if (i*i+j*j+k*k > (R-N_smoothing)*(R-N_smoothing))
    			A3D_ELEM(pMask, k, i, j) = -1;
    	}
    	int Rparticle = round(sqrt(radius));
    	std::cout << "particle radius = " << Rparticle << std::endl;
    	size_t xrows = angles.mdimx;

    	resolutionMatrix.initConstant(xrows, NVoxelsOriginalMask, maxRes);


    	#ifdef DEBUG_MASK
    	std::cout << "-------------DEBUG-----------" <<std::endl;
    	std::cout << "Next number ought to be the same than number of directions"
    			<< std::endl;
    	std::cout << "number of angles" << xrows << std::endl;
    	std::cout << "---------END-DEBUG--------" <<std::endl;
    	mask.write("mask.vol");
    	#endif

    	maxRes = 18;
    	minRes = 1;
    	V.clear();

    	double u;

		freq_fourier_z.initZeros(ZSIZE(fftV));
		freq_fourier_x.initZeros(XSIZE(fftV));
		freq_fourier_y.initZeros(YSIZE(fftV));

		VEC_ELEM(freq_fourier_z,0) = 1e-38;
		for(size_t k=0; k<ZSIZE(fftV); ++k)
		{
			FFT_IDX2DIGFREQ(k,ZSIZE(pMask), u);
			VEC_ELEM(freq_fourier_z,k) = u;
		}

		VEC_ELEM(freq_fourier_y,0) = 1e-38;
		for(size_t k=0; k<YSIZE(fftV); ++k)
		{
			FFT_IDX2DIGFREQ(k,YSIZE(pMask), u);
			VEC_ELEM(freq_fourier_y,k) = u;
		}

		VEC_ELEM(freq_fourier_x,0) = 1e-38;
		for(size_t k=0; k<XSIZE(fftV); ++k)
		{
			FFT_IDX2DIGFREQ(k,XSIZE(pMask), u);
			VEC_ELEM(freq_fourier_x,k) = u;
		}
}


void ProgDirSharpening::directionalResolutionStep(Matrix2D<int> &faces, Matrix2D<double> vertex)
{
	std::cout << "Analyzing directions " << std::endl;

}

void ProgDirSharpening::run()
{
	//Defining general information to be used
	produceSideInfo();

	//Defining the number of vertex and faces of the icosahedron
	Matrix2D<double> vertex;
	Matrix2D<int> faces;

	icosahedronVertex(vertex);
	icosahedronFaces(faces, vertex);
	defineIcosahedronFaceMask(faces, vertex, fftV);
}


void ProgDirSharpening::defineIcosahedronFaceMask(Matrix2D<int> &faces, Matrix2D<double> &vertex,
		MultidimArray< std::complex<double> > &myfftV)
{

	double ang_con = 40.0*PI/180;

	for (size_t f = 0; f<MAT_YSIZE(faces); ++f)
	{
		//Repeated faces are skipped
		if (MAT_ELEM(faces, f, 0) == -1)
			continue;

		int v1, v2, v3;
		double x1, x2, x3, y1, y2, y3, z1, z2, z3;
		v1 = MAT_ELEM(faces,f, 0); v2 = MAT_ELEM(faces,f, 0); v3 = MAT_ELEM(faces,f, 0);

		x1 = MAT_ELEM(vertex,v1, 0); y1 = MAT_ELEM(vertex,v1, 1); z1 = MAT_ELEM(vertex,v1, 2);
		x2 = MAT_ELEM(vertex,v2, 0); y2 = MAT_ELEM(vertex,v2, 1); z2 = MAT_ELEM(vertex,v2, 2);
		x3 = MAT_ELEM(vertex,v2, 0); y3 = MAT_ELEM(vertex,v3, 1); z3 = MAT_ELEM(vertex,v3, 2);

		std::cout << x1 << " " << y1 << " " << z1 << std::endl;
		std::cout << x2 << " " << y2 << " " << z2 << std::endl;
		std::cout << x3 << " " << y3 << " " << z3 << std::endl;

		MultidimArray<int> conefilter;
		conefilter.resizeNoCopy(myfftV);

		double uz, uy, ux, uz_v1, uy_v1, ux_v1, uz_v2, uy_v2, ux_v2, uz_v3, uy_v3, ux_v3;
		long n = 0;
		for(size_t k=0; k<ZSIZE(myfftV); ++k)
		{
			uz = VEC_ELEM(freq_fourier_z, k);

			uz_v1 = uz*z1;
			uz_v2 = uz*z2;
			uz_v3 = uz*z3;

			for(size_t i=0; i<YSIZE(myfftV); ++i)
			{
				uy = VEC_ELEM(freq_fourier_y, i);

				uy_v1 = uz*y1;
				uy_v2 = uz*y2;
				uy_v3 = uz*y3;

				for(size_t j=0; j<XSIZE(myfftV); ++j)
				{

					double iun=DIRECT_MULTIDIM_ELEM(iu,n);
					ux = VEC_ELEM(freq_fourier_x, j);

					ux_v1 = ux*x1;
					ux_v2 = ux*x2;
					ux_v3 = ux*x3;

					//BE CAREFULL with the order
					//double dotproduct = (uy*x_dir + ux*y_dir + uz*z_dir)*iun;

					double acosine_v1 = acos(fabs((ux_v1 + uy_v1 + uz_v1)*iun));
					double acosine_v2 = acos(fabs((ux_v2 + uy_v2 + uz_v2)*iun));
					double acosine_v3 = acos(fabs((ux_v3 + uy_v3 + uz_v3)*iun));

					if ((acosine_v1<ang_con) && (acosine_v2<ang_con) && (acosine_v3<ang_con))
					{
						DIRECT_MULTIDIM_ELEM(conefilter, n) = 1;
					}
					else{
						DIRECT_MULTIDIM_ELEM(conefilter, n) = 0;
					}
					++n;
				}
			}


		}


		std::cout << "Analyzing directions " << std::endl;

		Image<int> icosahedronMasked;
		icosahedronMasked = conefilter;
		FileName fnmasked;
		fnmasked = formatString("masked_%i.mrc",f);
		icosahedronMasked.write(fnmasked);
	}


	exit(0);

}
