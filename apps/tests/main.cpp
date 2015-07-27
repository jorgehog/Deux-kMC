#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

#include "testdiffusion.h"

namespace {

// Tests that the uniform RNG does what it is supposed to.
TEST(MiscTests, uRNG)
{
    double mean = 0;
    double squareMean = 0;
    uint N = 10000000;

    for (uint i = 0; i < N; ++i)
    {
        double urng = rng.uniform();

        EXPECT_GE(1, urng);
        EXPECT_LE(0, urng);

        mean += urng;
        squareMean += urng*urng;
    }

    mean /= N;
    squareMean /= N;

    double var = squareMean - mean*mean;

    EXPECT_NEAR(1./12, var, 1E-2);
    EXPECT_NEAR(0.5, mean, 1E-3);
}

// Tests that the normal RNG does what it is supposed to.
TEST(MiscTests, nRNG)
{
    double mean = 0;
    double squareMean = 0;
    uint N = 10000000;

    for (uint i = 0; i < N; ++i)
    {
        double urng = rng.normal();

        mean += urng;
        squareMean += urng*urng;
    }

    mean /= N;
    squareMean /= N;

    double var = squareMean - mean*mean;

    EXPECT_NEAR(1, var, 1E-3);
    EXPECT_NEAR(0, mean, 1E-3);
}


}  // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    rng.initialize(time(NULL));

    return RUN_ALL_TESTS();
}
