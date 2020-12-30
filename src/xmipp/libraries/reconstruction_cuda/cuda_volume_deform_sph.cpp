// Xmipp includes
#include "api/dimension_vector.h"
#include "core/metadata_label.h"
#include "core/xmipp_random_mode.h"
#include "enum/argument_access_type.h"
#include "enum/argument_memory_location.h"
#include "enum/compute_api.h"
#include "ktt_types.h"
#include "reconstruction_adapt_cuda/volume_deform_sph_gpu.h"
#include "cuda_volume_deform_sph.h"
#include "core/matrix1d.h"
// Standard includes
#include <iterator>
#include <stdexcept>
#include <stdio.h>
#include <iostream>
#include <exception>
// Thrust includes
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
// KTT includes
#include "tuner_api.h"
// Cuda kernel include
//#include "cuda_volume_deform_sph.cu"

// CUDA kernel defines
#define BLOCK_X_DIM 8
#define BLOCK_Y_DIM 4
#define BLOCK_Z_DIM 4
#define TOTAL_BLOCK_SIZE (BLOCK_X_DIM * BLOCK_Y_DIM * BLOCK_Z_DIM)


// Not everything will be needed to transfer every time.
// Some parameters stay the same for the whole time.
//----------------------------------------------------
// constant params:
//      Rmax2, iRmax, VI, VR, vL1, vN, vL2, vM,
//      volumesI, volumesR
//----------------------------------------------------
// changing params:
//      steps, clnm, 
//----------------------------------------------------
// parameters that can be initialized at gpu:
//      outputs(diff2,sumVD,modg,Ncount) = 0
//      VO().initZeros(VR()).setXmippOrigin()
//      Gx().initZeros(VR()).setXmippOrigin(), Gy..., Gz...
//----------------------------------------------------
// applyTransformation is true only in the very last call.
// saveDeformation is true only in the very last call and only when analyzeStrain is true

// explicit instantiations
//template class VolumeDeformSph<float>;
//template class VolumeDeformSph<ComputationDataType>;

// Common functions
template<typename T>
cudaError cudaMallocAndCopy(T** target, const T* source, size_t numberOfElements, size_t memSize = 0) 
{
    size_t elemSize = numberOfElements * sizeof(T);
    memSize = memSize == 0 ? elemSize : memSize * sizeof(T);

    cudaError err = cudaSuccess;
    if ((err = cudaMalloc(target, memSize)) != cudaSuccess) {
        *target = NULL;
        return err;
    }

    if ((err = cudaMemcpy(*target, source, elemSize, cudaMemcpyHostToDevice)) != cudaSuccess) {
        cudaFree(*target);
        *target = NULL;
    }

    if (memSize > elemSize) {
        cudaMemset((*target) + numberOfElements, 0, memSize - elemSize);
    }

    return err;
}

void printCudaError() 
{
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess)
        fprintf(stderr, "Cuda error: %s\n", cudaGetErrorString(err));
}

// Copies data from CPU to the GPU and at the same time transforms from
// type 'U' to type 'T'. Works only for numeric types
template<typename Target, typename Source>
void transformData(Target** dest, Source* source, size_t n, bool mallocMem = true)
{
    std::vector<Target> tmp(source, source + n);

    if (mallocMem){
        if (cudaMalloc(dest, sizeof(Target) * n) != cudaSuccess){
            printCudaError();
        }
    }

    if (cudaMemcpy(*dest, tmp.data(), sizeof(Target) * n, cudaMemcpyHostToDevice) != cudaSuccess){
        printCudaError();
    }
}
/*
// Copies data from CPU to the GPU and at the same time transforms from
// type 'U' to type 'T'. Works only for numeric types
template<typename Target, typename Source>
void transformData(Target** dest, Source* source, size_t n, bool mallocMem = true)
{
    size_t aligned = n + COPY_BLOCK_X_DIM - (n % COPY_BLOCK_X_DIM);

    Source* gpuSource;
    if (cudaMallocAndCopy(&gpuSource, source, n, aligned) != cudaSuccess) {
        printCudaError();
    }

    if (mallocMem){
        if (cudaMalloc(dest, aligned * sizeof(Target)) != cudaSuccess) {
            printCudaError();
        }
    }

    transformAndCopyKernel<<<aligned / COPY_BLOCK_X_DIM, COPY_BLOCK_X_DIM>>>(*dest, gpuSource);
    cudaDeviceSynchronize();

    cudaFree(gpuSource);
}
*/
// VolumeDeformSph methods

VolumeDeformSph::VolumeDeformSph() : tuner(0, 0, ktt::ComputeAPI::CUDA)
{
}

VolumeDeformSph::~VolumeDeformSph() 
{
    freeImage(images.VI);
    freeImage(images.VR);
    freeImage(images.VO);

    cudaFree(zshparams.vL1);
    cudaFree(zshparams.vL2);
    cudaFree(zshparams.vN);
    cudaFree(zshparams.vM);

    for (size_t i = 0; i < volumes.size; i++) {
        freeImage(justForFreeR[i]);
        freeImage(justForFreeI[i]);
    }
    cudaFree(volumes.R);
    cudaFree(volumes.I);

    cudaFree(steps);
    cudaFree(clnm);

    cudaFree(outputs);

    freeImage(deformImages.Gx);
    freeImage(deformImages.Gy);
    freeImage(deformImages.Gz);
}

void VolumeDeformSph::freeImage(ImageData &im) 
{
    if (im.data != nullptr)
        cudaFree(im.data);
}

void VolumeDeformSph::associateWith(ProgVolumeDeformSphGpu* prog) 
{
    program = prog;
}

void VolumeDeformSph::setupConstantParameters() 
{
    if (program == nullptr)
        throw new std::runtime_error("VolumeDeformSph not associated with the program!");

    this->Rmax2 = program->Rmax * program->Rmax;
    this->iRmax = 1 / program->Rmax;
    setupImage(program->VI, images.VI);
    setupImage(program->VR, images.VR);
    setupZSHparams();
    setupVolumes();

    // ktt stuff
    Rmax2Id = tuner.addArgumentScalar(Rmax2);
    iRmaxId = tuner.addArgumentScalar(iRmax);
    imagesId = tuner.addArgumentScalar(images);
    zshparamsId = tuner.addArgumentScalar(zshparams);
    volumesId = tuner.addArgumentScalar(volumes);
    // this one is just a dummy argument, for real it is initilized
    // in the setupChangingParameters, but it has to be initilized here
    // (or elsewhere) for successful kernel call
    deformImagesId = tuner.addArgumentScalar(deformImages);

    // kernel dimension
    kttBlock.setSizeX(BLOCK_X_DIM);
    kttBlock.setSizeY(BLOCK_Y_DIM);
    kttBlock.setSizeZ(BLOCK_Z_DIM);
    kttGrid.setSizeX(images.VR.xDim / BLOCK_X_DIM);
    kttGrid.setSizeY(images.VR.yDim / BLOCK_Y_DIM);
    kttGrid.setSizeZ(images.VR.zDim / BLOCK_Z_DIM);

    // kernel init
    kernelId = tuner.addKernelFromFile(pathToXmipp + pathToKernel, "computeDeform", kttGrid, kttBlock);
}

void VolumeDeformSph::setupChangingParameters() 
{
    if (program == nullptr)
        throw new std::runtime_error("VolumeDeformSph not associated with the program!");

    unsigned stepsSize = program->steps_cp.size() * sizeof(ComputationDataType);
    unsigned clnmSize = program->clnm.size() * sizeof(ComputationDataType);

    if (this->steps == nullptr) {
        if (cudaMalloc(&(this->steps), stepsSize) != cudaSuccess)
            printCudaError();
        else
            // ktt stuff
            stepsId = tuner.addArgumentVector<ComputationDataType>(static_cast<ktt::UserBuffer>(steps), zshparams.size * sizeof(ComputationDataType), ktt::ArgumentAccessType::ReadOnly, ktt::ArgumentMemoryLocation::Device);
    }
    if (this->clnm == nullptr) {
        if (cudaMalloc(&(this->clnm), clnmSize) != cudaSuccess)
            printCudaError();
        else
            // ktt stuff
            clnmId = tuner.addArgumentVector<ComputationDataType>(static_cast<ktt::UserBuffer>(clnm), zshparams.size * sizeof(ComputationDataType), ktt::ArgumentAccessType::ReadOnly, ktt::ArgumentMemoryLocation::Device);
    }

    transformData(&(this->steps), program->steps_cp.vdata, program->steps_cp.size(), false);
    transformData(&(this->clnm), program->clnm.vdata, program->clnm.size(), false);

    this->applyTransformation = program->applyTransformation;
    this->saveDeformation = program->saveDeformation;

    // ktt stuff
    applyTransformationId = tuner.addArgumentScalar(static_cast<int>(applyTransformation));
    saveDeformationId = tuner.addArgumentScalar(static_cast<int>(saveDeformation));

    if (applyTransformation) {
        setupImage(images.VR, images.VO);
        imagesId = tuner.addArgumentScalar(images);
    }
    if (saveDeformation) {
        setupImage(images.VR, deformImages.Gx);
        setupImage(images.VR, deformImages.Gy);
        setupImage(images.VR, deformImages.Gz);
        deformImagesId = tuner.addArgumentScalar(deformImages);
    }
}

KernelOutputs VolumeDeformSph::getOutputs() 
{
    return exOuts;
}

void VolumeDeformSph::transferImageData(Image<double>& outputImage, ImageData& inputData) 
{
    size_t elements = inputData.xDim * inputData.yDim * inputData.zDim;
    std::vector<ComputationDataType> tVec(elements);
    cudaMemcpy(tVec.data(), inputData.data, sizeof(ComputationDataType) * elements, cudaMemcpyDeviceToHost);
    std::vector<double> dVec(tVec.begin(), tVec.end());
    memcpy(outputImage().data, dVec.data(), sizeof(double) * elements);
    /*
    size_t size = inputData.xDim * inputData.yDim * inputData.zDim * sizeof(T);
    cudaMemcpy(outputImage().data, inputData.data, size, cudaMemcpyDeviceToHost);
    */
    /*
    double* tmp;
    transformData(&tmp, inputData.data, elements);
    cudaMemcpy(outputImage().data, tmp, elements * sizeof(double), cudaMemcpyDeviceToHost);
    */
}

void VolumeDeformSph::runKernel() 
{
    // Does not work in general case, but test data have nice sizes

// KTT test
    // Define path to kernel
    //std::string pathToXmipp = "/home/david/thesis/xmipp-bundle/";
    //std::string pathToKernel = "src/xmipp/libraries/reconstruction_cuda/cuda_volume_deform_sph.cu";

    // Define block and grid size
    //const ktt::DimensionVector kttBlock(BLOCK_X_DIM, BLOCK_Y_DIM, BLOCK_Z_DIM);
    //const ktt::DimensionVector kttGrid(images.VR.xDim / BLOCK_X_DIM, images.VR.yDim / BLOCK_Y_DIM, images.VR.zDim / BLOCK_Z_DIM);

    // Define thrust reduction vector
    thrust::device_vector<ComputationDataType> t_out(kttGrid.getTotalSize() * 4, 0.0);

    // Initialize tuner
    //tuner.setCompilerOptions("--std=c++11");

    // Add kernel to the tuner
    //ktt::KernelId kernelId = tuner.addKernelFromFile(pathToXmipp + pathToKernel, "computeDeform", kttGrid, kttBlock);

    // Add arguments for the kernel
    //Rmax2Id = tuner.addArgumentScalar(Rmax2);
    //iRmaxId = tuner.addArgumentScalar(iRmax);
    //imagesId = tuner.addArgumentScalar(images);
    //zshparamsId = tuner.addArgumentScalar(zshparams);
    //volumesId = tuner.addArgumentScalar(volumes);
    //deformImagesId = tuner.addArgumentScalar(deformImages);
    //stepsId = tuner.addArgumentVector<ComputationDataType>(static_cast<ktt::UserBuffer>(steps), zshparams.size * sizeof(ComputationDataType), ktt::ArgumentAccessType::ReadOnly, ktt::ArgumentMemoryLocation::Device);
    //clnmId = tuner.addArgumentVector<ComputationDataType>(static_cast<ktt::UserBuffer>(clnm), zshparams.size * sizeof(ComputationDataType), ktt::ArgumentAccessType::ReadOnly, ktt::ArgumentMemoryLocation::Device);
    // Bool is not supported, will see what happens
    //applyTransformationId = tuner.addArgumentScalar(static_cast<int>(applyTransformation));
    //saveDeformationId = tuner.addArgumentScalar(static_cast<int>(saveDeformation));
    // end of booleans
    ktt::ArgumentId thrustVecId = tuner.addArgumentVector<ComputationDataType>(static_cast<ktt::UserBuffer>(thrust::raw_pointer_cast(t_out.data())), t_out.size() * sizeof(ComputationDataType), ktt::ArgumentAccessType::ReadWrite, ktt::ArgumentMemoryLocation::Device);

    // Assign arguments to the kernel
    tuner.setKernelArguments(kernelId, std::vector<ktt::ArgumentId>{
            Rmax2Id, iRmaxId, imagesId, zshparamsId,
            stepsId, clnmId,
            volumesId, deformImagesId, applyTransformationId, saveDeformationId,
            thrustVecId
            });

    // Run kernel
    tuner.runKernel(kernelId, {}, {});
// end KTT test

    cudaDeviceSynchronize();

    auto diff2It = t_out.begin();
    auto sumVDIt = diff2It + kttGrid.getTotalSize();
    auto modgIt = sumVDIt + kttGrid.getTotalSize();
    auto NcountIt = modgIt + kttGrid.getTotalSize();

    exOuts.diff2 = thrust::reduce(diff2It, sumVDIt);
    exOuts.sumVD = thrust::reduce(sumVDIt, modgIt);
    exOuts.modg = thrust::reduce(modgIt, NcountIt);
    exOuts.Ncount = thrust::reduce(NcountIt, t_out.end());
}

void VolumeDeformSph::transferResults() 
{
    if (applyTransformation) {
        transferImageData(program->VO, images.VO);
    }
    if (saveDeformation) {
        transferImageData(program->Gx, deformImages.Gx);
        transferImageData(program->Gy, deformImages.Gy);
        transferImageData(program->Gz, deformImages.Gz);
    }
}

void VolumeDeformSph::setupZSHparams()
{
    zshparams.size = program->vL1.size();

    if (cudaMallocAndCopy(&zshparams.vL1, program->vL1.vdata, zshparams.size) != cudaSuccess)
        printCudaError();
    if (cudaMallocAndCopy(&zshparams.vL2, program->vL2.vdata, zshparams.size) != cudaSuccess)
        printCudaError();
    if (cudaMallocAndCopy(&zshparams.vN, program->vN.vdata, zshparams.size) != cudaSuccess)
        printCudaError();
    if (cudaMallocAndCopy(&zshparams.vM, program->vM.vdata, zshparams.size) != cudaSuccess)
        printCudaError();
}

void VolumeDeformSph::setupVolumes()
{
    volumes.size = program->volumesR.size();

    justForFreeR.resize(volumes.size);
    justForFreeI.resize(volumes.size);

    for (size_t i = 0; i < volumes.size; i++) {
        setupImage(program->volumesR[i], justForFreeR[i]);
        setupImage(program->volumesI[i], justForFreeI[i]);
    }

    if (cudaMallocAndCopy(&volumes.R, justForFreeR.data(), volumes.size) != cudaSuccess)
        printCudaError();
    if (cudaMallocAndCopy(&volumes.I, justForFreeI.data(), volumes.size) != cudaSuccess)
        printCudaError();
}

void VolumeDeformSph::setupImage(Image<double>& inputImage, ImageData& outputImageData) 
{
    auto& mda = inputImage();

    outputImageData.xShift = mda.xinit;
    outputImageData.yShift = mda.yinit;
    outputImageData.zShift = mda.zinit;
    outputImageData.xDim = mda.xdim;
    outputImageData.yDim = mda.ydim;
    outputImageData.zDim = mda.zdim;

    // if T is smaller than double -> error
    // might be replaced with cudaMallocAndCopy, but there are different sizes: T vs double
    /*
    int size = outputImageData.xDim * outputImageData.yDim * outputImageData.zDim * sizeof(T);
    if (cudaMalloc(&outputImageData.data, size) != cudaSuccess)
        printCudaError();
    if (cudaMemcpy(outputImageData.data, mda.data, size, cudaMemcpyHostToDevice) != cudaSuccess)
        printCudaError();
    */
    transformData(&outputImageData.data, mda.data, mda.xdim * mda.ydim * mda.zdim);
}

void VolumeDeformSph::setupImage(ImageData& inputImage, ImageData& outputImageData, bool copyData) 
{
    outputImageData.xShift = inputImage.xShift;
    outputImageData.yShift = inputImage.yShift;
    outputImageData.zShift = inputImage.zShift;
    outputImageData.xDim = inputImage.xDim;
    outputImageData.yDim = inputImage.yDim;
    outputImageData.zDim = inputImage.zDim;

    size_t size = inputImage.xDim * inputImage.yDim * inputImage.zDim * sizeof(ComputationDataType);
    if (cudaMalloc(&outputImageData.data, size) != cudaSuccess)
        printCudaError();

    if (copyData) {
        if (cudaMemcpy(outputImageData.data, inputImage.data, size, cudaMemcpyHostToDevice) != cudaSuccess)
            printCudaError();
    }
}
