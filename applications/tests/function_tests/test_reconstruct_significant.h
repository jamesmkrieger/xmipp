//#include <gtest/gtest.h>

#include "reconstruction/reconstruct_significant.h"

class RecSignTest : public ProgReconstructSignificant , public ::testing::Test {
public:

    void SetUp() {
        //get example images/staks
        if (chdir(((String)(getXmippPath() + (String)"/resources/test/recSignTest")).c_str())==-1)
            REPORT_ERROR(ERR_UNCLASSIFIED,"Could not change directory");
    }


    void alignImagesToGallery() {
        // nastavit mdIn na images.xmd
        size_t Nvols = this->mdGallery.size();
        size_t Ndirs = this->mdGallery[0].size();
        this->cc.initZeros(Nimgs,Nvols,Ndirs);

        this->fnIn =
    }
};
