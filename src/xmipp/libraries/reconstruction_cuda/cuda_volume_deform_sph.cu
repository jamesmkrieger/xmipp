#ifndef CUDA_VOLUME_DEFORM_SPH_CU
#define CUDA_VOLUME_DEFORM_SPH_CU
//#include "cuda_volume_deform_sph.h"

// Compilation settings

#ifdef USE_DOUBLE_PRECISION
// Types
using PrecisionType = double;
using PrecisionType3 = double3;
// Constants
#define _PI_ (3.1415926535897931e+0)
// Functions
#define SQRT sqrt
#define ATAN2 atan2
#define COS cos
#define SIN sin

#else
// Types
using PrecisionType = float;
using PrecisionType3 = float3;
// Constants
#define _PI_ (3.1415926535897f)
// Functions
#define SQRT sqrtf
#define ATAN2 atan2f
#define COS cosf
#define SIN sinf

#endif// USE_DOUBLE_PRECISION

#ifdef USE_SCATTERED_ZSH_CLNM
using ClnmType = PrecisionType*;
using ZshParamsType =
    struct ZSHparams { int *vL1, *vN, *vL2, *vM; unsigned size; };
#else
using ClnmType = PrecisionType3*;
using ZshParamsType = int4*;
#endif// USE_SCATTERED_ZSH_CLNM

// Compilation settings - end

// Define used data structures

struct ImageData
{
    int xShift;
    int yShift;
    int zShift;

    int xDim;
    int yDim;
    int zDim;

    PrecisionType* data;
};
/*
#ifdef USE_SCATTERED_ZSH_CLNM
struct ZSHparams 
{
    int* vL1;
    int* vN;
    int* vL2;
    int* vM;
    unsigned size;
};
#endif
*/
struct Volumes 
{
    ImageData* I;
    ImageData* R;
    unsigned size;
};

struct IROimages 
{
    ImageData VI;
    ImageData VR;
    ImageData VO;
};

struct DeformImages 
{
    ImageData Gx;
    ImageData Gy;
    ImageData Gz;
};

// CUDA kernel defines
#define BLOCK_X_DIM 8
#define BLOCK_Y_DIM 4
#define BLOCK_Z_DIM 4
#define TOTAL_BLOCK_SIZE (BLOCK_X_DIM * BLOCK_Y_DIM * BLOCK_Z_DIM)

// ImageData macros

// Index to global memory
#define GET_IDX(ImD,k,i,j) \
    ((ImD).xDim * (ImD).yDim * (k) + (ImD).xDim * (i) + (j))

// Logical index = Physical index + shift
#define P2L_X_IDX(ImD,j) \
    ((j) + (ImD).xShift)

#define P2L_Y_IDX(ImD,i) \
    ((i) + (ImD).yShift)

#define P2L_Z_IDX(ImD,k) \
    ((k) + (ImD).zShift)

// Physical index = Logical index - shift
#define L2P_X_IDX(ImD,j) \
    ((j) - (ImD).xShift)

#define L2P_Y_IDX(ImD,i) \
    ((i) - (ImD).yShift)

#define L2P_Z_IDX(ImD,k) \
    ((k) - (ImD).zShift)

// Element access
#define ELEM_3D(ImD,k,i,j) \
    ((ImD).data[GET_IDX((ImD), (k), (i), (j))])

#define ELEM_3D_SHIFTED(ImD,k,i,j) \
    (ELEM_3D((ImD), (k) - (ImD).zShift, (i) - (ImD).yShift, (j) - (ImD).xShift))

// Utility macros
#define MY_OUTSIDE(ImD,k,i,j) \
    ((j) < (ImD).xShift || (j) > (ImD).xShift + (ImD).xDim - 1 || \
     (i) < (ImD).yShift || (i) > (ImD).yShift + (ImD).yDim - 1 || \
     (k) < (ImD).zShift || (k) > (ImD).zShift + (ImD).zDim - 1)

// Smart casting to selected precision (at compile time)
// ...just shorter static_cast
#define CST(num) (static_cast<PrecisionType>((num)))

#define FLOOR(x) (((x) == (int)(x)) ? (int)(x):(((x) > 0) ? (int)(x) : \
                  (int)((x) - 1)))
#define LIN_INTERP(a, l, h) ((l) + ((h) - (l)) * (a))

// Forward declarations
__device__ PrecisionType ZernikeSphericalHarmonics(int l1, int n, int l2, int m,
        PrecisionType xr, PrecisionType yr, PrecisionType zr, PrecisionType r);

__device__ PrecisionType interpolatedElement3D(ImageData ImD,
        PrecisionType x, PrecisionType y, PrecisionType z,
        PrecisionType doutside_value = 0);

/*
 * The beast
 */
extern "C" __global__ void computeDeform(
        PrecisionType Rmax2,
        PrecisionType iRmax,
        IROimages images,
        ZshParamsType zshparams,
        ClnmType clnm,
        int steps,
        Volumes volumes,
        DeformImages deformImages,
        bool applyTransformation,
        bool saveDeformation,
        PrecisionType* g_outArr
        ) 
{
    __shared__ PrecisionType sumArray[TOTAL_BLOCK_SIZE * 4];

    // Compute thread index in a block
    int tIdx = threadIdx.z * blockDim.x * blockDim.y + threadIdx.y * blockDim.x + threadIdx.x;

    // Get physical indexes
    int kPhys = blockIdx.z * blockDim.z + threadIdx.z;
    int iPhys = blockIdx.y * blockDim.y + threadIdx.y;
    int jPhys = blockIdx.x * blockDim.x + threadIdx.x;

    // Update to logical indexes (calculations expect logical indexing)
    int k = P2L_Z_IDX(images.VR, kPhys);
    int i = P2L_Y_IDX(images.VR, iPhys);
    int j = P2L_X_IDX(images.VR, jPhys);

    PrecisionType r2 = k*k + i*i + j*j;
    PrecisionType rr = SQRT(r2) * iRmax;
    PrecisionType gx = 0.0, gy = 0.0, gz = 0.0;

    if (r2 < Rmax2) {
        for (int idx = 0; idx < steps; idx++) {
#ifdef USE_SCATTERED_ZSH_CLNM
            int l1 = zshparams.vL1[idx];
            int n = zshparams.vN[idx];
            int l2 = zshparams.vL2[idx];
            int m = zshparams.vM[idx];
#else
            int l1 = zshparams[idx].w;
            int n = zshparams[idx].x;
            int l2 = zshparams[idx].y;
            int m = zshparams[idx].z;
#endif
#ifdef USE_ZSH_FUNCTION
            PrecisionType zsph = ZernikeSphericalHarmonics(l1, n, l2, m,
                    j * iRmax, i * iRmax, k * iRmax, rr);
#else
            PrecisionType xr = j * iRmax, yr = i * iRmax, zr = k * iRmax;

            // General variables
            PrecisionType r2 = rr * rr, xr2 = xr * xr, yr2 = yr * yr,
                          zr2 = zr * zr;

#if L2 >= 5
            // Variables needed for l2 >= 5
            PrecisionType tht = CST(0.0), phi = CST(0.0), cost = CST(0.0),
                          sint = CST(0.0), cost2 = CST(0.0), sint2 = CST(0.0);
            if (l2 >= 5) {
              tht = ATAN2(yr, xr);
              phi = ATAN2(zr, SQRT(xr2 + yr2));
              sint = SIN(phi);
              cost = COS(tht);
              sint2 = sint * sint;
              cost2 = cost * cost;
            }
#endif// L2 >= 5

            // Zernike polynomial
            PrecisionType R = CST(0.0);

            switch (l1) {
            case 0:
              R = SQRT(CST(3));
              break;
            case 1:
              R = SQRT(CST(5)) * rr;
              break;
            case 2:
              switch (n) {
              case 0:
                R = CST(-0.5) * SQRT(CST(7)) *
                    (CST(2.5) * (1 - 2 * r2) + CST(0.5));
                break;
              case 2:
                R = SQRT(CST(7)) * r2;
                break;
              }
              break;
#if L1 >= 3
            case 3:
              switch (n) {
              case 1:
                R = CST(-1.5) * rr * (CST(3.5) * (1 - 2 * r2) + CST(1.5));
                break;
              case 3:
                R = 3 * r2 * rr;
              }
              break;
#endif// L1 >= 3
#if L1 >= 4
            case 4:
              switch (n) {
              case 0:
                R = SQRT(CST(11)) *
                    ((63 * r2 * r2 / 8) - (35 * r2 / 4) + (CST(15) / CST(8)));
                break;
              case 2:
                R = CST(-0.5) * SQRT(CST(11)) * r2 *
                    (CST(4.5) * (1 - 2 * r2) + CST(2.5));
                break;
              case 4:
                R = SQRT(CST(11)) * r2 * r2;
                break;
              }
              break;
#endif// L1 >= 4
#if L1 >= 5
            case 5:
              switch (n) {
              case 1:
                R = SQRT(CST(13)) * rr *
                    ((99 * r2 * r2 / 8) - (63 * r2 / 4) + (CST(35) / CST(8)));
                break;
              case 3:
                R = CST(-0.5) * SQRT(CST(13)) * r2 * rr *
                    (CST(5.5) * (1 - 2 * r2) + CST(3.5));
                break;
              }
              break;
#endif// L1 >= 5
            }

            // Spherical harmonic
            PrecisionType Y = CST(0.0);

            switch (l2) {
            case 0:
              Y = (CST(1.0) / CST(2.0)) * SQRT((PrecisionType)CST(1.0) / _PI_);
              break;
            case 1:
              switch (m) {
              case -1:
                Y = SQRT(CST(3.0) / (CST(4.0) * _PI_)) * yr;
                break;
              case 0:
                Y = SQRT(CST(3.0) / (CST(4.0) * _PI_)) * zr;
                break;
              case 1:
                Y = SQRT(CST(3.0) / (CST(4.0) * _PI_)) * xr;
                break;
              }
              break;
            case 2:
              switch (m) {
              case -2:
                Y = SQRT(CST(15.0) / (CST(4.0) * _PI_)) * xr * yr;
                break;
              case -1:
                Y = SQRT(CST(15.0) / (CST(4.0) * _PI_)) * zr * yr;
                break;
              case 0:
                Y = SQRT(CST(5.0) / (CST(16.0) * _PI_)) *
                    (-xr2 - yr2 + CST(2.0) * zr2);
                break;
              case 1:
                Y = SQRT(CST(15.0) / (CST(4.0) * _PI_)) * xr * zr;
                break;
              case 2:
                Y = SQRT(CST(15.0) / (CST(16.0) * _PI_)) * (xr2 - yr2);
                break;
              }
              break;
#if L2 >= 3
            case 3:
              switch (m) {
              case -3:
                Y = SQRT(CST(35.0) / (CST(16.0) * CST(2.0) * _PI_)) * yr *
                    (CST(3.0) * xr2 - yr2);
                break;
              case -2:
                Y = SQRT(CST(105.0) / (CST(4.0) * _PI_)) * zr * yr * xr;
                break;
              case -1:
                Y = SQRT(CST(21.0) / (CST(16.0) * CST(2.0) * _PI_)) * yr *
                    (CST(4.0) * zr2 - xr2 - yr2);
                break;
              case 0:
                Y = SQRT(CST(7.0) / (CST(16.0) * _PI_)) * zr *
                    (CST(2.0) * zr2 - CST(3.0) * xr2 - CST(3.0) * yr2);
                break;
              case 1:
                Y = SQRT(CST(21.0) / (CST(16.0) * CST(2.0) * _PI_)) * xr *
                    (CST(4.0) * zr2 - xr2 - yr2);
                break;
              case 2:
                Y = SQRT(CST(105.0) / (CST(16.0) * _PI_)) * zr * (xr2 - yr2);
                break;
              case 3:
                Y = SQRT(CST(35.0) / (CST(16.0) * CST(2.0) * _PI_)) * xr *
                    (xr2 - CST(3.0) * yr2);
                break;
              }
              break;
#endif// L2 >= 3
#if L2 >= 4
            case 4:
              switch (m) {
              case -4:
                Y = SQRT((CST(35.0) * CST(9.0)) / (CST(16.0) * _PI_)) * yr *
                    xr * (xr2 - yr2);
                break;
              case -3:
                Y = SQRT((CST(9.0) * CST(35.0)) /
                         (CST(16.0) * CST(2.0) * _PI_)) *
                    yr * zr * (CST(3.0) * xr2 - yr2);
                break;
              case -2:
                Y = SQRT((CST(9.0) * CST(5.0)) / (CST(16.0) * _PI_)) * yr * xr *
                    (CST(7.0) * zr2 - (xr2 + yr2 + zr2));
                break;
              case -1:
                Y = SQRT((CST(9.0) * CST(5.0)) /
                         (CST(16.0) * CST(2.0) * _PI_)) *
                    yr * zr * (CST(7.0) * zr2 - CST(3.0) * (xr2 + yr2 + zr2));
                break;
              case 0:
                Y = SQRT(CST(9.0) / (CST(16.0) * CST(16.0) * _PI_)) *
                    (CST(35.0) * zr2 * zr2 - CST(30.0) * zr2 + CST(3.0));
                break;
              case 1:
                Y = SQRT((CST(9.0) * CST(5.0)) /
                         (CST(16.0) * CST(2.0) * _PI_)) *
                    xr * zr * (CST(7.0) * zr2 - CST(3.0) * (xr2 + yr2 + zr2));
                break;
              case 2:
                Y = SQRT((CST(9.0) * CST(5.0)) / (CST(8.0) * CST(8.0) * _PI_)) *
                    (xr2 - yr2) * (CST(7.0) * zr2 - (xr2 + yr2 + zr2));
                break;
              case 3:
                Y = SQRT((CST(9.0) * CST(35.0)) /
                         (CST(16.0) * CST(2.0) * _PI_)) *
                    xr * zr * (xr2 - CST(3.0) * yr2);
                break;
              case 4:
                Y = SQRT((CST(9.0) * CST(35.0)) /
                         (CST(16.0) * CST(16.0) * _PI_)) *
                    (xr2 * (xr2 - CST(3.0) * yr2) -
                     yr2 * (CST(3.0) * xr2 - yr2));
                break;
              }
              break;
#endif// L2 >= 4
#if L2 >= 5
            case 5:
              switch (m) {
              case -5:
                Y = (CST(3.0) / CST(16.0)) *
                    SQRT(CST(77.0) / (CST(2.0) * _PI_)) * sint2 * sint2 * sint *
                    SIN(CST(5.0) * phi);
                break;
              case -4:
                Y = (CST(3.0) / CST(8.0)) *
                    SQRT(CST(385.0) / (CST(2.0) * _PI_)) * sint2 * sint2 *
                    SIN(CST(4.0) * phi);
                break;
              case -3:
                Y = (CST(1.0) / CST(16.0)) *
                    SQRT(CST(385.0) / (CST(2.0) * _PI_)) * sint2 * sint *
                    (CST(9.0) * cost2 - CST(1.0)) * SIN(CST(3.0) * phi);
                break;
              case -2:
                Y = (CST(1.0) / CST(4.0)) *
                    SQRT(CST(1155.0) / (CST(4.0) * _PI_)) * sint2 *
                    (CST(3.0) * cost2 * cost - cost) * SIN(CST(2.0) * phi);
                break;
              case -1:
                Y = (CST(1.0) / CST(8.0)) *
                    SQRT(CST(165.0) / (CST(4.0) * _PI_)) * sint *
                    (CST(21.0) * cost2 * cost2 - CST(14.0) * cost2 + 1) *
                    SIN(phi);
                break;
              case 0:
                Y = (CST(1.0) / CST(16.0)) * SQRT(CST(11.0) / _PI_) *
                    (CST(63.0) * cost2 * cost2 * cost -
                     CST(70.0) * cost2 * cost + CST(15.0) * cost);
                break;
              case 1:
                Y = (CST(1.0) / CST(8.0)) *
                    SQRT(CST(165.0) / (CST(4.0) * _PI_)) * sint *
                    (CST(21.0) * cost2 * cost2 - CST(14.0) * cost2 + 1) *
                    COS(phi);
                break;
              case 2:
                Y = (CST(1.0) / CST(4.0)) *
                    SQRT(CST(1155.0) / (CST(4.0) * _PI_)) * sint2 *
                    (CST(3.0) * cost2 * cost - cost) * COS(CST(2.0) * phi);
                break;
              case 3:
                Y = (CST(1.0) / CST(16.0)) *
                    SQRT(CST(385.0) / (CST(2.0) * _PI_)) * sint2 * sint *
                    (CST(9.0) * cost2 - CST(1.0)) * COS(CST(3.0) * phi);
                break;
              case 4:
                Y = (CST(3.0) / CST(8.0)) *
                    SQRT(CST(385.0) / (CST(2.0) * _PI_)) * sint2 * sint2 *
                    COS(CST(4.0) * phi);
                break;
              case 5:
                Y = (CST(3.0) / CST(16.0)) *
                    SQRT(CST(77.0) / (CST(2.0) * _PI_)) * sint2 * sint2 * sint *
                    COS(CST(5.0) * phi);
                break;
              }
              break;
#endif// L2 >= 5
            }

            PrecisionType zsph = R * Y;
#endif// USE_ZSH_FUNCTION

            if (rr > 0 || l2 == 0) {
#ifdef USE_SCATTERED_ZSH_CLNM
                gx += zsph * clnm[idx];
                gy += zsph * clnm[idx + zshparams.size];
                gz += zsph * clnm[idx + zshparams.size * 2];
#else
                gx += zsph * clnm[idx].x;
                gy += zsph * clnm[idx].y;
                gz += zsph * clnm[idx].z;
#endif
            }
        }
    }

    PrecisionType voxelI, voxelR;
    PrecisionType diff;

    PrecisionType localDiff2 = 0.0, localSumVD = 0.0, localModg = 0.0, localNcount = 0.0;

    if (applyTransformation) {
        // Indexing requires physical indexes
        voxelR = ELEM_3D(images.VR, kPhys, iPhys, jPhys);
        // Logical indexes used to check whether the point is in the matrix
        voxelI = interpolatedElement3D(images.VI, j + gx, i + gy, k + gz);

        if (voxelI >= 0.0)
            localSumVD += voxelI;

        ELEM_3D(images.VO, kPhys, iPhys, jPhys) = voxelI;
        diff = voxelR - voxelI;
        localDiff2 += diff * diff;
        localModg += gx*gx + gy*gy + gz*gz;
        localNcount++;
    }

    for (unsigned idv = 0; idv < volumes.size; idv++) {
        voxelR = ELEM_3D(volumes.R[idv], kPhys, iPhys, jPhys);
        voxelI = interpolatedElement3D(volumes.I[idv], j + gx, i + gy, k + gz);

        if (voxelI >= 0.0)
            localSumVD += voxelI;

        diff = voxelR - voxelI;
        localDiff2 += diff * diff;
        localModg += gx*gx + gy*gy + gz*gz;
        localNcount++;
    }

    sumArray[tIdx] = localDiff2;
    sumArray[tIdx + TOTAL_BLOCK_SIZE] = localSumVD;
    sumArray[tIdx + TOTAL_BLOCK_SIZE * 2] = localModg;
    sumArray[tIdx + TOTAL_BLOCK_SIZE * 3] = localNcount;

    __syncthreads();
    
    // Block reduction   
    for (int s = TOTAL_BLOCK_SIZE / 2; s > 0; s /= 2) {
        if (tIdx < s) {
            sumArray[tIdx] += sumArray[tIdx + s];
            sumArray[tIdx + TOTAL_BLOCK_SIZE] += sumArray[tIdx + TOTAL_BLOCK_SIZE + s];
            sumArray[tIdx + TOTAL_BLOCK_SIZE * 2] += sumArray[tIdx + TOTAL_BLOCK_SIZE * 2 + s];
            sumArray[tIdx + TOTAL_BLOCK_SIZE * 3] += sumArray[tIdx + TOTAL_BLOCK_SIZE * 3 + s];
        }
        __syncthreads();
    }

    // Save values to the global memory for later
    if (tIdx == 0) {
        int bIdx = blockIdx.z * gridDim.x * gridDim.y + blockIdx.y * gridDim.x + blockIdx.x;
        int TOTAL_GRID_SIZE = gridDim.x * gridDim.y * gridDim.z;
        g_outArr[bIdx] = sumArray[0];
        g_outArr[bIdx + TOTAL_GRID_SIZE] = sumArray[TOTAL_BLOCK_SIZE];
        g_outArr[bIdx + TOTAL_GRID_SIZE * 2] = sumArray[TOTAL_BLOCK_SIZE * 2];
        g_outArr[bIdx + TOTAL_GRID_SIZE * 3] = sumArray[TOTAL_BLOCK_SIZE * 3];
    }

    if (saveDeformation) {
        ELEM_3D(deformImages.Gx, kPhys, iPhys, jPhys) = gx;
        ELEM_3D(deformImages.Gy, kPhys, iPhys, jPhys) = gy;
        ELEM_3D(deformImages.Gz, kPhys, iPhys, jPhys) = gz;
    }
}

/*
 * Linear interpolation
 */
__device__ PrecisionType interpolatedElement3D(ImageData ImD,
        PrecisionType x, PrecisionType y, PrecisionType z,
        PrecisionType outside_value) 
{
        int x0 = FLOOR(x);
        PrecisionType fx = x - x0;
        int x1 = x0 + 1;

        int y0 = FLOOR(y);
        PrecisionType fy = y - y0;
        int y1 = y0 + 1;

        int z0 = FLOOR(z);
        PrecisionType fz = z - z0;
        int z1 = z0 + 1;

        PrecisionType d000 = (MY_OUTSIDE(ImD, z0, y0, x0)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z0, y0, x0);
        PrecisionType d001 = (MY_OUTSIDE(ImD, z0, y0, x1)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z0, y0, x1);
        PrecisionType d010 = (MY_OUTSIDE(ImD, z0, y1, x0)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z0, y1, x0);
        PrecisionType d011 = (MY_OUTSIDE(ImD, z0, y1, x1)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z0, y1, x1);
        PrecisionType d100 = (MY_OUTSIDE(ImD, z1, y0, x0)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z1, y0, x0);
        PrecisionType d101 = (MY_OUTSIDE(ImD, z1, y0, x1)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z1, y0, x1);
        PrecisionType d110 = (MY_OUTSIDE(ImD, z1, y1, x0)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z1, y1, x0);
        PrecisionType d111 = (MY_OUTSIDE(ImD, z1, y1, x1)) ?
            outside_value : ELEM_3D_SHIFTED(ImD, z1, y1, x1);

        PrecisionType dx00 = LIN_INTERP(fx, d000, d001);
        PrecisionType dx01 = LIN_INTERP(fx, d100, d101);
        PrecisionType dx10 = LIN_INTERP(fx, d010, d011);
        PrecisionType dx11 = LIN_INTERP(fx, d110, d111);
        PrecisionType dxy0 = LIN_INTERP(fy, dx00, dx10);
        PrecisionType dxy1 = LIN_INTERP(fy, dx01, dx11);

        return LIN_INTERP(fz, dxy0, dxy1);
}

/*
 * ZSH
 */
__device__ PrecisionType ZernikeSphericalHarmonics(int l1, int n, int l2, int m, PrecisionType xr, PrecisionType yr, PrecisionType zr, PrecisionType r)
{
	// General variables
	PrecisionType r2=r*r,xr2=xr*xr,yr2=yr*yr,zr2=zr*zr;

	//Variables needed for l>=5
	PrecisionType tht=CST(0.0),phi=CST(0.0),cost=CST(0.0),sint=CST(0.0),cost2=CST(0.0),sint2=CST(0.0);
	if (l2>=5)
	{
		tht = ATAN2(yr,xr);
		phi = ATAN2(zr,SQRT(xr2 + yr2));
		sint = SIN(phi); cost = COS(tht);
		sint2 = sint*sint; cost2 = cost*cost;
	}

	// Zernike polynomial
	PrecisionType R=CST(0.0);

	switch (l1)
	{
	case 0:
		R = SQRT(CST(3));
		break;
	case 1:
		R = SQRT(CST(5))*r;
		break;
	case 2:
		switch (n)
		{
		case 0:
			R = CST(-0.5)*SQRT(CST(7))*(CST(2.5)*(1-2*r2)+CST(0.5));
			break;
		case 2:
			R = SQRT(CST(7))*r2;
			break;
		} break;
	case 3:
		switch (n)
		{
		case 1:
			R = CST(-1.5)*r*(CST(3.5)*(1-2*r2)+CST(1.5));
			break;
		case 3:
			R = 3*r2*r;
		} break;
	case 4:
		switch (n)
		{
		case 0:
			R = SQRT(CST(11))*((63*r2*r2/8)-(35*r2/4)+(CST(15)/CST(8)));
			break;
		case 2:
			R = CST(-0.5)*SQRT(CST(11))*r2*(CST(4.5)*(1-2*r2)+CST(2.5));
			break;
		case 4:
			R = SQRT(CST(11))*r2*r2;
			break;
		} break;
	case 5:
		switch (n)
		{
		case 1:
			R = SQRT(CST(13))*r*((99*r2*r2/8)-(63*r2/4)+(CST(35)/CST(8)));
			break;
		case 3:
			R = CST(-0.5)*SQRT(CST(13))*r2*r*(CST(5.5)*(1-2*r2)+CST(3.5));
			break;
		} break;
	}

	// Spherical harmonic
	PrecisionType Y=CST(0.0);

	switch (l2)
	{
	case 0:
		Y = (CST(1.0)/CST(2.0))*SQRT((PrecisionType) CST(1.0)/_PI_);
		break;
	case 1:
		switch (m)
		{
		case -1:
			Y = SQRT(CST(3.0)/(CST(4.0)*_PI_))*yr;
			break;
		case 0:
			Y = SQRT(CST(3.0)/(CST(4.0)*_PI_))*zr;
			break;
		case 1:
			Y = SQRT(CST(3.0)/(CST(4.0)*_PI_))*xr;
			break;
		} break;
	case 2:
		switch (m)
		{
		case -2:
			Y = SQRT(CST(15.0)/(CST(4.0)*_PI_))*xr*yr;
			break;
		case -1:
			Y = SQRT(CST(15.0)/(CST(4.0)*_PI_))*zr*yr;
			break;
		case 0:
			Y = SQRT(CST(5.0)/(CST(16.0)*_PI_))*(-xr2-yr2+CST(2.0)*zr2);
			break;
		case 1:
			Y = SQRT(CST(15.0)/(CST(4.0)*_PI_))*xr*zr;
			break;
		case 2:
			Y = SQRT(CST(15.0)/(CST(16.0)*_PI_))*(xr2-yr2);
			break;
		} break;
	case 3:
		switch (m)
		{
		case -3:
			Y = SQRT(CST(35.0)/(CST(16.0)*CST(2.0)*_PI_))*yr*(CST(3.0)*xr2-yr2);
			break;
		case -2:
			Y = SQRT(CST(105.0)/(CST(4.0)*_PI_))*zr*yr*xr;
			break;
		case -1:
			Y = SQRT(CST(21.0)/(CST(16.0)*CST(2.0)*_PI_))*yr*(CST(4.0)*zr2-xr2-yr2);
			break;
		case 0:
			Y = SQRT(CST(7.0)/(CST(16.0)*_PI_))*zr*(CST(2.0)*zr2-CST(3.0)*xr2-CST(3.0)*yr2);
			break;
		case 1:
			Y = SQRT(CST(21.0)/(CST(16.0)*CST(2.0)*_PI_))*xr*(CST(4.0)*zr2-xr2-yr2);
			break;
		case 2:
			Y = SQRT(CST(105.0)/(CST(16.0)*_PI_))*zr*(xr2-yr2);
			break;
		case 3:
			Y = SQRT(CST(35.0)/(CST(16.0)*CST(2.0)*_PI_))*xr*(xr2-CST(3.0)*yr2);
			break;
		} break;
	case 4:
		switch (m)
		{
		case -4:
			Y = SQRT((CST(35.0)*CST(9.0))/(CST(16.0)*_PI_))*yr*xr*(xr2-yr2);
			break;
		case -3:
			Y = SQRT((CST(9.0)*CST(35.0))/(CST(16.0)*CST(2.0)*_PI_))*yr*zr*(CST(3.0)*xr2-yr2);
			break;
		case -2:
			Y = SQRT((CST(9.0)*CST(5.0))/(CST(16.0)*_PI_))*yr*xr*(CST(7.0)*zr2-(xr2+yr2+zr2));
			break;
		case -1:
			Y = SQRT((CST(9.0)*CST(5.0))/(CST(16.0)*CST(2.0)*_PI_))*yr*zr*(CST(7.0)*zr2-CST(3.0)*(xr2+yr2+zr2));
			break;
		case 0:
			Y = SQRT(CST(9.0)/(CST(16.0)*CST(16.0)*_PI_))*(CST(35.0)*zr2*zr2-CST(30.0)*zr2+CST(3.0));
			break;
		case 1:
			Y = SQRT((CST(9.0)*CST(5.0))/(CST(16.0)*CST(2.0)*_PI_))*xr*zr*(CST(7.0)*zr2-CST(3.0)*(xr2+yr2+zr2));
			break;
		case 2:
			Y = SQRT((CST(9.0)*CST(5.0))/(CST(8.0)*CST(8.0)*_PI_))*(xr2-yr2)*(CST(7.0)*zr2-(xr2+yr2+zr2));
			break;
		case 3:
			Y = SQRT((CST(9.0)*CST(35.0))/(CST(16.0)*CST(2.0)*_PI_))*xr*zr*(xr2-CST(3.0)*yr2);
			break;
		case 4:
			Y = SQRT((CST(9.0)*CST(35.0))/(CST(16.0)*CST(16.0)*_PI_))*(xr2*(xr2-CST(3.0)*yr2)-yr2*(CST(3.0)*xr2-yr2));
			break;
		} break;
	case 5:
		switch (m)
		{
		case -5:
			Y = (CST(3.0)/CST(16.0))*SQRT(CST(77.0)/(CST(2.0)*_PI_))*sint2*sint2*sint*SIN(CST(5.0)*phi);
			break;
		case -4:
			Y = (CST(3.0)/CST(8.0))*SQRT(CST(385.0)/(CST(2.0)*_PI_))*sint2*sint2*SIN(CST(4.0)*phi);
			break;
		case -3:
			Y = (CST(1.0)/CST(16.0))*SQRT(CST(385.0)/(CST(2.0)*_PI_))*sint2*sint*(CST(9.0)*cost2-CST(1.0))*SIN(CST(3.0)*phi);
			break;
		case -2:
			Y = (CST(1.0)/CST(4.0))*SQRT(CST(1155.0)/(CST(4.0)*_PI_))*sint2*(CST(3.0)*cost2*cost-cost)*SIN(CST(2.0)*phi);
			break;
		case -1:
			Y = (CST(1.0)/CST(8.0))*SQRT(CST(165.0)/(CST(4.0)*_PI_))*sint*(CST(21.0)*cost2*cost2-CST(14.0)*cost2+1)*SIN(phi);
			break;
		case 0:
			Y = (CST(1.0)/CST(16.0))*SQRT(CST(11.0)/_PI_)*(CST(63.0)*cost2*cost2*cost-CST(70.0)*cost2*cost+CST(15.0)*cost);
			break;
		case 1:
			Y = (CST(1.0)/CST(8.0))*SQRT(CST(165.0)/(CST(4.0)*_PI_))*sint*(CST(21.0)*cost2*cost2-CST(14.0)*cost2+1)*COS(phi);
			break;
		case 2:
			Y = (CST(1.0)/CST(4.0))*SQRT(CST(1155.0)/(CST(4.0)*_PI_))*sint2*(CST(3.0)*cost2*cost-cost)*COS(CST(2.0)*phi);
			break;
		case 3:
			Y = (CST(1.0)/CST(16.0))*SQRT(CST(385.0)/(CST(2.0)*_PI_))*sint2*sint*(CST(9.0)*cost2-CST(1.0))*COS(CST(3.0)*phi);
			break;
		case 4:
			Y = (CST(3.0)/CST(8.0))*SQRT(CST(385.0)/(CST(2.0)*_PI_))*sint2*sint2*COS(CST(4.0)*phi);
			break;
		case 5:
			Y = (CST(3.0)/CST(16.0))*SQRT(CST(77.0)/(CST(2.0)*_PI_))*sint2*sint2*sint*COS(CST(5.0)*phi);
			break;
		}break;
	}

	return R*Y;
}


// Function redefinition
#ifdef USE_DOUBLE_PRECISION 
/*
__device__ double ZernikeSphericalHarmonics(int l1, int n, int l2, int m, double xr, double yr, double zr, double r)
{
	// General variables
	double r2=r*r,xr2=xr*xr,yr2=yr*yr,zr2=zr*zr;

	//Variables needed for l>=5
	double tht=0.0,phi=0.0,cost=0.0,sint=0.0,cost2=0.0,sint2=0.0;
	if (l2>=5)
	{
		tht = atan2(yr,xr);
		phi = atan2(zr,sqrt(xr2 + yr2));
		sint = sin(phi); cost = cos(tht);
		sint2 = sint*sint; cost2 = cost*cost;
	}

	// Zernike polynomial
	double R=0.0;

	switch (l1)
	{
	case 0:
		R = sqrt((double) 3);
		break;
	case 1:
		R = sqrt((double) 5)*r;
		break;
	case 2:
		switch (n)
		{
		case 0:
			R = -0.5*sqrt((double) 7)*(2.5*(1-2*r2)+0.5);
			break;
		case 2:
			R = sqrt((double) 7)*r2;
			break;
		} break;
	case 3:
		switch (n)
		{
		case 1:
			R = -1.5*r*(3.5*(1-2*r2)+1.5);
			break;
		case 3:
			R = 3*r2*r;
		} break;
	case 4:
		switch (n)
		{
		case 0:
			R = sqrt((double) 11)*((63*r2*r2/8)-(35*r2/4)+(15/8));
			break;
		case 2:
			R = -0.5*sqrt((double) 11)*r2*(4.5*(1-2*r2)+2.5);
			break;
		case 4:
			R = sqrt((double) 11)*r2*r2;
			break;
		} break;
	case 5:
		switch (n)
		{
		case 1:
			R = sqrt((double) 13)*r*((99*r2*r2/8)-(63*r2/4)+(35/8));
			break;
		case 3:
			R = -0.5*sqrt((double) 13)*r2*r*(5.5*(1-2*r2)+3.5);
			break;
		} break;
	}

	// Spherical harmonic
	double Y=0.0;

	switch (l2)
	{
	case 0:
		Y = (1.0/2.0)*sqrt((double) 1.0/_PI_);
		break;
	case 1:
		switch (m)
		{
		case -1:
			Y = sqrt(3.0/(4.0*_PI_))*yr;
			break;
		case 0:
			Y = sqrt(3.0/(4.0*_PI_))*zr;
			break;
		case 1:
			Y = sqrt(3.0/(4.0*_PI_))*xr;
			break;
		} break;
	case 2:
		switch (m)
		{
		case -2:
			Y = sqrt(15.0/(4.0*_PI_))*xr*yr;
			break;
		case -1:
			Y = sqrt(15.0/(4.0*_PI_))*zr*yr;
			break;
		case 0:
			Y = sqrt(5.0/(16.0*_PI_))*(-xr2-yr2+2.0*zr2);
			break;
		case 1:
			Y = sqrt(15.0/(4.0*_PI_))*xr*zr;
			break;
		case 2:
			Y = sqrt(15.0/(16.0*_PI_))*(xr2-yr2);
			break;
		} break;
	case 3:
		switch (m)
		{
		case -3:
			Y = sqrt(35.0/(16.0*2.0*_PI_))*yr*(3.0*xr2-yr2);
			break;
		case -2:
			Y = sqrt(105.0/(4.0*_PI_))*zr*yr*xr;
			break;
		case -1:
			Y = sqrt(21.0/(16.0*2.0*_PI_))*yr*(4.0*zr2-xr2-yr2);
			break;
		case 0:
			Y = sqrt(7.0/(16.0*_PI_))*zr*(2.0*zr2-3.0*xr2-3.0*yr2);
			break;
		case 1:
			Y = sqrt(21.0/(16.0*2.0*_PI_))*xr*(4.0*zr2-xr2-yr2);
			break;
		case 2:
			Y = sqrt(105.0/(16.0*_PI_))*zr*(xr2-yr2);
			break;
		case 3:
			Y = sqrt(35.0/(16.0*2.0*_PI_))*xr*(xr2-3.0*yr2);
			break;
		} break;
	case 4:
		switch (m)
		{
		case -4:
			Y = sqrt((35.0*9.0)/(16.0*_PI_))*yr*xr*(xr2-yr2);
			break;
		case -3:
			Y = sqrt((9.0*35.0)/(16.0*2.0*_PI_))*yr*zr*(3.0*xr2-yr2);
			break;
		case -2:
			Y = sqrt((9.0*5.0)/(16.0*_PI_))*yr*xr*(7.0*zr2-(xr2+yr2+zr2));
			break;
		case -1:
			Y = sqrt((9.0*5.0)/(16.0*2.0*_PI_))*yr*zr*(7.0*zr2-3.0*(xr2+yr2+zr2));
			break;
		case 0:
			Y = sqrt(9.0/(16.0*16.0*_PI_))*(35.0*zr2*zr2-30.0*zr2+3.0);
			break;
		case 1:
			Y = sqrt((9.0*5.0)/(16.0*2.0*_PI_))*xr*zr*(7.0*zr2-3.0*(xr2+yr2+zr2));
			break;
		case 2:
			Y = sqrt((9.0*5.0)/(8.0*8.0*_PI_))*(xr2-yr2)*(7.0*zr2-(xr2+yr2+zr2));
			break;
		case 3:
			Y = sqrt((9.0*35.0)/(16.0*2.0*_PI_))*xr*zr*(xr2-3.0*yr2);
			break;
		case 4:
			Y = sqrt((9.0*35.0)/(16.0*16.0*_PI_))*(xr2*(xr2-3.0*yr2)-yr2*(3.0*xr2-yr2));
			break;
		} break;
	case 5:
		switch (m)
		{
		case -5:
			Y = (3.0/16.0)*sqrt(77.0/(2.0*_PI_))*sint2*sint2*sint*sin(5.0*phi);
			break;
		case -4:
			Y = (3.0/8.0)*sqrt(385.0/(2.0*_PI_))*sint2*sint2*sin(4.0*phi);
			break;
		case -3:
			Y = (1.0/16.0)*sqrt(385.0/(2.0*_PI_))*sint2*sint*(9.0*cost2-1.0)*sin(3.0*phi);
			break;
		case -2:
			Y = (1.0/4.0)*sqrt(1155.0/(4.0*_PI_))*sint2*(3.0*cost2*cost-cost)*sin(2.0*phi);
			break;
		case -1:
			Y = (1.0/8.0)*sqrt(165.0/(4.0*_PI_))*sint*(21.0*cost2*cost2-14.0*cost2+1)*sin(phi);
			break;
		case 0:
			Y = (1.0/16.0)*sqrt(11.0/_PI_)*(63.0*cost2*cost2*cost-70.0*cost2*cost+15.0*cost);
			break;
		case 1:
			Y = (1.0/8.0)*sqrt(165.0/(4.0*_PI_))*sint*(21.0*cost2*cost2-14.0*cost2+1)*cos(phi);
			break;
		case 2:
			Y = (1.0/4.0)*sqrt(1155.0/(4.0*_PI_))*sint2*(3.0*cost2*cost-cost)*cos(2.0*phi);
			break;
		case 3:
			Y = (1.0/16.0)*sqrt(385.0/(2.0*_PI_))*sint2*sint*(9.0*cost2-1.0)*cos(3.0*phi);
			break;
		case 4:
			Y = (3.0/8.0)*sqrt(385.0/(2.0*_PI_))*sint2*sint2*cos(4.0*phi);
			break;
		case 5:
			Y = (3.0/16.0)*sqrt(77.0/(2.0*_PI_))*sint2*sint2*sint*cos(5.0*phi);
			break;
		}break;
	}

	return R*Y;
}
*/
#else
/*
__device__ float ZernikeSphericalHarmonics(int l1, int n, int l2, int m, float xr, float yr, float zr, float r)
{
	// General variables
	float r2=r*r,xr2=xr*xr,yr2=yr*yr,zr2=zr*zr;

	//Variables needed for l>=5
	float tht=0.0f,phi=0.0f,cost=0.0f,sint=0.0f,cost2=0.0f,sint2=0.0f;
	if (l2>=5)
	{
		tht = atan2f(yr,xr);
		phi = atan2f(zr,sqrtf(xr2 + yr2));
		sint = sinf(phi); cost = cosf(tht);
		sint2 = sint*sint; cost2 = cost*cost;
	}

	// Zernike polynomial
	float R=0.0f;

	switch (l1)
	{
	case 0:
		R = sqrtf((float) 3);
		break;
	case 1:
		R = sqrtf((float) 5)*r;
		break;
	case 2:
		switch (n)
		{
		case 0:
			R = -0.5f*sqrtf((float) 7)*(2.5f*(1-2*r2)+0.5f);
			break;
		case 2:
			R = sqrtf((float) 7)*r2;
			break;
		} break;
	case 3:
		switch (n)
		{
		case 1:
			R = -1.5f*r*(3.5f*(1-2*r2)+1.5f);
			break;
		case 3:
			R = 3*r2*r;
		} break;
	case 4:
		switch (n)
		{
		case 0:
			R = sqrtf((float) 11)*((63*r2*r2/8)-(35*r2/4)+(15/8));
			break;
		case 2:
			R = -0.5f*sqrtf((float) 11)*r2*(4.5f*(1-2*r2)+2.5f);
			break;
		case 4:
			R = sqrtf((float) 11)*r2*r2;
			break;
		} break;
	case 5:
		switch (n)
		{
		case 1:
			R = sqrtf((float) 13)*r*((99*r2*r2/8)-(63*r2/4)+(35/8));
			break;
		case 3:
			R = -0.5f*sqrtf((float) 13)*r2*r*(5.5f*(1-2*r2)+3.5f);
			break;
		} break;
	}

	// Spherical harmonic
	float Y=0.0f;

	switch (l2)
	{
	case 0:
		Y = (1.0f/2.0f)*sqrtf((float) 1.0f/_PI_);
		break;
	case 1:
		switch (m)
		{
		case -1:
			Y = sqrtf(3.0f/(4.0f*_PI_))*yr;
			break;
		case 0:
			Y = sqrtf(3.0f/(4.0f*_PI_))*zr;
			break;
		case 1:
			Y = sqrtf(3.0f/(4.0f*_PI_))*xr;
			break;
		} break;
	case 2:
		switch (m)
		{
		case -2:
			Y = sqrtf(15.0f/(4.0f*_PI_))*xr*yr;
			break;
		case -1:
			Y = sqrtf(15.0f/(4.0f*_PI_))*zr*yr;
			break;
		case 0:
			Y = sqrtf(5.0f/(16.0f*_PI_))*(-xr2-yr2+2.0f*zr2);
			break;
		case 1:
			Y = sqrtf(15.0f/(4.0f*_PI_))*xr*zr;
			break;
		case 2:
			Y = sqrtf(15.0f/(16.0f*_PI_))*(xr2-yr2);
			break;
		} break;
	case 3:
		switch (m)
		{
		case -3:
			Y = sqrtf(35.0f/(16.0f*2.0f*_PI_))*yr*(3.0f*xr2-yr2);
			break;
		case -2:
			Y = sqrtf(105.0f/(4.0f*_PI_))*zr*yr*xr;
			break;
		case -1:
			Y = sqrtf(21.0f/(16.0f*2.0f*_PI_))*yr*(4.0f*zr2-xr2-yr2);
			break;
		case 0:
			Y = sqrtf(7.0f/(16.0f*_PI_))*zr*(2.0f*zr2-3.0f*xr2-3.0f*yr2);
			break;
		case 1:
			Y = sqrtf(21.0f/(16.0f*2.0f*_PI_))*xr*(4.0f*zr2-xr2-yr2);
			break;
		case 2:
			Y = sqrtf(105.0f/(16.0f*_PI_))*zr*(xr2-yr2);
			break;
		case 3:
			Y = sqrtf(35.0f/(16.0f*2.0f*_PI_))*xr*(xr2-3.0f*yr2);
			break;
		} break;
	case 4:
		switch (m)
		{
		case -4:
			Y = sqrtf((35.0f*9.0f)/(16.0f*_PI_))*yr*xr*(xr2-yr2);
			break;
		case -3:
			Y = sqrtf((9.0f*35.0f)/(16.0f*2.0f*_PI_))*yr*zr*(3.0f*xr2-yr2);
			break;
		case -2:
			Y = sqrtf((9.0f*5.0f)/(16.0f*_PI_))*yr*xr*(7.0f*zr2-(xr2+yr2+zr2));
			break;
		case -1:
			Y = sqrtf((9.0f*5.0f)/(16.0f*2.0f*_PI_))*yr*zr*(7.0f*zr2-3.0f*(xr2+yr2+zr2));
			break;
		case 0:
			Y = sqrtf(9.0f/(16.0f*16.0f*_PI_))*(35.0f*zr2*zr2-30.0f*zr2+3.0f);
			break;
		case 1:
			Y = sqrtf((9.0f*5.0f)/(16.0f*2.0f*_PI_))*xr*zr*(7.0f*zr2-3.0f*(xr2+yr2+zr2));
			break;
		case 2:
			Y = sqrtf((9.0f*5.0f)/(8.0f*8.0f*_PI_))*(xr2-yr2)*(7.0f*zr2-(xr2+yr2+zr2));
			break;
		case 3:
			Y = sqrtf((9.0f*35.0f)/(16.0f*2.0f*_PI_))*xr*zr*(xr2-3.0f*yr2);
			break;
		case 4:
			Y = sqrtf((9.0f*35.0f)/(16.0f*16.0f*_PI_))*(xr2*(xr2-3.0f*yr2)-yr2*(3.0f*xr2-yr2));
			break;
		} break;
	case 5:
		switch (m)
		{
		case -5:
			Y = (3.0f/16.0f)*sqrtf(77.0f/(2.0f*_PI_))*sint2*sint2*sint*sinf(5.0f*phi);
			break;
		case -4:
			Y = (3.0f/8.0f)*sqrtf(385.0f/(2.0f*_PI_))*sint2*sint2*sinf(4.0f*phi);
			break;
		case -3:
			Y = (1.0f/16.0f)*sqrtf(385.0f/(2.0f*_PI_))*sint2*sint*(9.0f*cost2-1.0f)*sinf(3.0f*phi);
			break;
		case -2:
			Y = (1.0f/4.0f)*sqrtf(1155.0f/(4.0f*_PI_))*sint2*(3.0f*cost2*cost-cost)*sinf(2.0f*phi);
			break;
		case -1:
			Y = (1.0f/8.0f)*sqrtf(165.0f/(4.0f*_PI_))*sint*(21.0f*cost2*cost2-14.0f*cost2+1)*sinf(phi);
			break;
		case 0:
			Y = (1.0f/16.0f)*sqrtf(11.0f/_PI_)*(63.0f*cost2*cost2*cost-70.0f*cost2*cost+15.0f*cost);
			break;
		case 1:
			Y = (1.0f/8.0f)*sqrtf(165.0f/(4.0f*_PI_))*sint*(21.0f*cost2*cost2-14.0f*cost2+1)*cosf(phi);
			break;
		case 2:
			Y = (1.0f/4.0f)*sqrtf(1155.0f/(4.0f*_PI_))*sint2*(3.0f*cost2*cost-cost)*cosf(2.0f*phi);
			break;
		case 3:
			Y = (1.0f/16.0f)*sqrtf(385.0f/(2.0f*_PI_))*sint2*sint*(9.0f*cost2-1.0f)*cosf(3.0f*phi);
			break;
		case 4:
			Y = (3.0f/8.0f)*sqrtf(385.0f/(2.0f*_PI_))*sint2*sint2*cosf(4.0f*phi);
			break;
		case 5:
			Y = (3.0f/16.0f)*sqrtf(77.0f/(2.0f*_PI_))*sint2*sint2*sint*cosf(5.0f*phi);
			break;
		}break;
	}

	return R*Y;
}
*/
#endif//COMP_DOUBLE


#endif //CUDA_VOLUME_DEFORM_SPH_CU
