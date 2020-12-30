#ifndef VOLUME_DEFORM_SPH_H
#define VOLUME_DEFORM_SPH_H

#include <vector>
#include "api/dimension_vector.h"
#include "core/xmipp_image.h"
#include "core/multidim_array.h"

#include "ktt_types.h"
#include "tuner_api.h"

#ifdef COMP_DOUBLE
using ComputationDataType = double;
#else
using ComputationDataType = float;
#endif

// Forward declarations
class ProgVolumeDeformSphGpu;

struct ImageData
{
    int xShift = 0;
    int yShift = 0;
    int zShift = 0;

    int xDim = 0;
    int yDim = 0;
    int zDim = 0;

    ComputationDataType* data = nullptr;
};

struct ZSHparams 
{
    int* vL1 = nullptr;
    int* vN = nullptr;
    int* vL2 = nullptr;
    int* vM = nullptr;
    unsigned size = 0;
};

struct Volumes 
{
    ImageData* I = nullptr;
    ImageData* R = nullptr;
    unsigned size = 0;
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

struct KernelOutputs 
{
    ComputationDataType diff2 = 0.0;
    ComputationDataType sumVD = 0.0;
    ComputationDataType modg = 0.0;
    ComputationDataType Ncount = 0.0;
};

class VolumeDeformSph
{
public:
    void associateWith(ProgVolumeDeformSphGpu* prog);
    void setupConstantParameters();
    void setupChangingParameters();

    void runKernel();
    void transferResults();

    KernelOutputs getOutputs();
    void transferImageData(Image<double>& outputImage, ImageData& inputData);

    VolumeDeformSph();
    ~VolumeDeformSph();

private:
    ProgVolumeDeformSphGpu* program = nullptr;

    // ktt stuff
    ktt::Tuner tuner;
    ktt::KernelId kernelId;

    ktt::DimensionVector kttBlock;
    ktt::DimensionVector kttGrid;

    // Kernel path
    const std::string pathToXmipp = "/home/david/thesis/xmipp-bundle/";
    const std::string pathToKernel = "src/xmipp/libraries/reconstruction_cuda/cuda_volume_deform_sph.cu";

    // Variables transfered to the GPU memory
    ktt::ArgumentId Rmax2Id;
    ComputationDataType Rmax2;

    ktt::ArgumentId iRmaxId;
    ComputationDataType iRmax;

    ktt::ArgumentId stepsId;
    ComputationDataType* steps = nullptr;

    ktt::ArgumentId clnmId;
    ComputationDataType* clnm = nullptr;

    ktt::ArgumentId applyTransformationId;
    bool applyTransformation;

    ktt::ArgumentId saveDeformationId;
    bool saveDeformation;

    // Inside pointers point to the GPU memory

    ktt::ArgumentId imagesId;
    IROimages images;

    ktt::ArgumentId deformImagesId;
    DeformImages deformImages;

    ktt::ArgumentId zshparamsId;
    ZSHparams zshparams;

    ktt::ArgumentId volumesId;
    Volumes volumes;
    // because of the stupid design... :(
    std::vector<ImageData> justForFreeI;
    std::vector<ImageData> justForFreeR;

    KernelOutputs *outputs = nullptr;
    KernelOutputs exOuts;

    // helper methods for simplifying and transfering data to gpu

    void setupImage(Image<double>& inputImage, ImageData& outputImageData);
    void setupImage(ImageData& inputImage, ImageData& outputImageData, bool copyData = false);

    void freeImage(ImageData &im);

    void simplifyVec(std::vector<Image<double>>& vec, std::vector<ImageData>& res);

    void setupVolumes();

    void setupZSHparams();
};

#endif// VOLUME_DEFORM_SPH_H
