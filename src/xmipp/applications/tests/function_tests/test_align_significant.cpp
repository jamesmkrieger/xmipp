#include <iostream>

#include <gtest/gtest.h>
#include "alignment_test_utils.h"
#include <data/dimensions.h>
#include <core/xmipp_image.h>
class AlignSignificantTests : public ::testing::Test
{
};

static FileName OUT_DIR = "testIn";

void fillRow(MDRow &row, std::string &image, double angle, int refIndex) {
    row.setValue(MDL_IMAGE, image);
    row.setValue(MDL_ENABLED, 1);
    row.setValue(MDL_ANGLE_ROT, (double)angle);
    row.setValue(MDL_ANGLE_TILT, (double)0);
    row.setValue(MDL_ANGLE_PSI, (double)0);
    row.setValue(MDL_SHIFT_X, (double)0); // store negative translation
    row.setValue(MDL_SHIFT_Y, (double)0); // store negative translation
    row.setValue(MDL_FLIP, false);
    row.setValue(MDL_REF, refIndex);
}

TEST_F( AlignSignificantTests, test1)
{XMIPP_TRY
    int step = 10;
    int count = 360 / step;
    auto d = Dimensions(64, 64, 1, count);
    auto imgs = Image<float>(d.x(), d.y(), d.z(), d.n());
    for (size_t i = 0; i < d.n(); ++i) {
        auto addr = imgs.data.data + i * d.sizeSingle();
        drawClockArms(addr, d, d.x() / 2, d.y() / 2, (float)i * step);
    }
    imgs.write(OUT_DIR + "/armClockRotatedImgs.stk");

    d = Dimensions(64, 64, 1, 2);
    auto ref = Image<float>(d.x(), d.y(), d.z(), d.n());
    auto md = MetaData();
    MDRow row;
    auto fnStk = FileName(OUT_DIR + "/armClockRotatedRef.stk");
    auto fnXmd = FileName(OUT_DIR + "/armClockRotatedRef.xmd");
    for (size_t i = 0; i < d.n(); ++i) {
        auto addr = ref.data.data + i * d.sizeSingle();
        fnXmd.compose(i + 1, fnStk);
        fillRow(row, fnXmd , 0, 0);
        md.addRow(row);
        if (0 == i) {
            drawClockArms(addr, d, d.x() / 2, d.y() / 2, (float)i * step);
        }
    }
    ref.write(fnStk);
    md.write(OUT_DIR + "/armClockRotatedRef.xmd");


XMIPP_CATCH
}

TEST_F( AlignSignificantTests, test2)
{XMIPP_TRY
    int step = 10;
    int count = 360 / step;
    auto d = Dimensions(64, 64, 1, count);
    auto clock = Image<float>(d.x(), d.y());
    drawClockArms(clock.data.data, d, d.x() / 2, d.y() / 2, 0);
    auto imgs = Image<float>(d.x(), d.y(), d.z(), d.n());
    for (size_t i = 0; i < d.n(); ++i) {
        auto sX = 0;
        if (i < 18) {
            sX += -9 + i;
        }
        auto sY = 0;
        if (i >= 18) {
            sY += -27 + i;
        }
        auto m = Matrix2D<float>();
        m.initIdentity(3);
        MAT_ELEM(m,0,2) += sX;
        MAT_ELEM(m,1,2) += sY;
        auto r = Matrix2D<float>();
        rotation2DMatrix((float)i*step, r);
        m = r * m;
        auto out = MultidimArray<float>(1, d.z(), d.y(), d.x(), imgs.data.data + (i * d.sizeSingle()));
        auto in = MultidimArray<float>(1, d.z(), d.y(), d.x(), clock.data.data);
        out.setXmippOrigin();
        in.setXmippOrigin();
        applyGeometry(1, out, in, m, false, DONT_WRAP);
    }
    imgs.write(OUT_DIR + "/armClockRSImgs.stk");

    d = Dimensions(64, 64, 1, 2);
    auto ref = Image<float>(d.x(), d.y(), d.z(), d.n());
    auto md = MetaData();
    MDRow row;
    auto fnStk = FileName(OUT_DIR + "/armClockRSRef.stk");
    auto fnXmd = FileName(OUT_DIR + "/armClockRSRef.xmd");
    for (size_t i = 0; i < d.n(); ++i) {
        auto addr = ref.data.data + i * d.sizeSingle();
        fnXmd.compose(i + 1, fnStk);
        fillRow(row, fnXmd , 0, 0);
        md.addRow(row);
        if (0 == i) {
            memcpy(addr, clock.data.data, d.sizeSingle() * sizeof(float));
        }
    }
    ref.write(fnStk);
    md.write(OUT_DIR + "/armClockRSRef.xmd");


XMIPP_CATCH
}


GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
