#include <gtest/gtest.h>
#include <random>
#include "reconstruction/arotation_estimator.h"
#include "core/transformations.h"
#include "core/xmipp_image.h"

template<typename T>
class ARotationEstimator_Test : public ::testing::Test {
public:
    SETUP

    void TearDown() {
        delete estimator;
    }

    SETUPTESTCASE

    static void TearDownTestCase() {
        delete hw;
    }

    void generateAndTest2D(size_t n, size_t batch) {
        std::uniform_int_distribution<> distSizeSmall(0, 368);
        std::uniform_int_distribution<> distSizeBig(369, 768);

        // only square inputs are valid:
        size_t size = distSizeSmall(mt);
        rotate2D(Dimensions(size, size, 1, n), batch);
        size = distSizeBig(mt);
        rotate2D(Dimensions(size, size, 1, n), batch);
    }

    void rotate2D(const Dimensions &dims, size_t batch)
    {
        using Alignment::AlignType;
        float maxRotation = 360.f - std::numeric_limits<float>::min();

        std::uniform_int_distribution<> distPos(0, dims.x());
        auto rotations = std::vector<float>();
        rotations.reserve(dims.n());
        std::uniform_real_distribution<> distRot(0, maxRotation);
        for(size_t n = 0; n < dims.n(); ++n) {
            rotations.emplace_back(distRot(mt));
        }

        auto others = new T[dims.size()]();
        auto ref = new T[dims.xy()]();
        T centerX = dims.x() / 2;
        T centerY = dims.y() / 2;
        drawClockArms(ref, dims, centerX, centerY, 0.f);

        for (size_t n = 0; n < dims.n(); ++n) {
            T *d = others + (n * dims.xyzPadded());
            drawClockArms(d, dims, centerX, centerY, rotations.at(n));
        }
//        outputData(others, dims);

        INIT
        estimator->load2DReferenceOneToN(ref);
        estimator->computeRotation2DOneToN(others);
        auto result = estimator->getRotations2D();

        EXPECT_EQ(rotations.size(), result.size());
        float maxError = 3.f; // degrees
        for (size_t n = 0; n < result.size(); ++n) {
            // we rotated by angle, so we should detect rotation in '360 - angle' degrees
            EXPECT_NEAR(360 - rotations.at(n), result.at(n), maxError);
        }

        delete[] others;
        delete[] ref;
    }

private:
    Alignment::ARotationEstimator<T> *estimator;
    static HW *hw;
    static std::mt19937 mt;

    void drawClockArms(T *result, const Dimensions &dims, size_t xPos, size_t yPos, float rotDegree) {
        size_t yArmSize = dims.y() - yPos;
        size_t xArmSize = dims.x() - xPos;

        MultidimArray<T> tmp(dims.y(), dims.x());
        for (size_t y = yPos; y < yPos + yArmSize; ++y) {
            size_t index = y * dims.x() + xPos;
            tmp.data[index] = 1;
        }

        for (size_t x = xPos; x < xPos + xArmSize; ++x) {
            size_t index = yPos * dims.x() + x;
            tmp.data[index] = 1;
        }
        MultidimArray<T> wrapper(1, 1, dims.y(), dims.x(), result);
        rotate(3, wrapper, tmp, rotDegree);
    }

    void outputData(T *data, const Dimensions &dims) {
        MultidimArray<T>wrapper(dims.n(), dims.z(), dims.y(), dims.x(), data);
        Image<T> img(wrapper);
        img.write("data.stk");
    }
};
TYPED_TEST_CASE_P(ARotationEstimator_Test);

template<typename T>
HW* ARotationEstimator_Test<T>::hw;
template<typename T>
std::mt19937 ARotationEstimator_Test<T>::mt(42); // fixed seed to ensure reproducibility


//***********************************************
//              Rotation tests
//***********************************************


TYPED_TEST_P( ARotationEstimator_Test, rotate2DOneToOne)
{
    // test one reference vs one image
    ARotationEstimator_Test<TypeParam>::generateAndTest2D(1, 1);
}

TYPED_TEST_P( ARotationEstimator_Test, rotate2DOneToMany)
{
    // check that n == batch works properly
    ARotationEstimator_Test<TypeParam>::generateAndTest2D(5, 1);
}

TYPED_TEST_P( ARotationEstimator_Test, rotate2DOneToManyBatched1)
{
    // test that n mod batch != 0 works
    ASSERT_THROW(ARotationEstimator_Test<TypeParam>::generateAndTest2D(5, 3), XmippError);
}

REGISTER_TYPED_TEST_CASE_P(ARotationEstimator_Test,
    rotate2DOneToOne,
    rotate2DOneToMany,
    rotate2DOneToManyBatched1
);
