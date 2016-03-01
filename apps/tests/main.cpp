#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

#include "testdiffusion.h"

#include "testboundaries.h"

#include "testcontinuumdiffusion.h"

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

TEST(MiscTests, BinarySearch)
{
    uint N = 1000;
    vec rates = ones(N);

    double R = N;
    vec cumsum = linspace(1, N, N);

    EXPECT_EQ(N, cumsum(N-1));

    uint nc = uint(1e7);

    double mean = 0;
    double mean2 = 0;

    for (uint n = 0; n < nc; ++n)
    {
        double c = chooseFromTotalRate(cumsum.memptr(), N, R)/double(N);

        mean += c;
        mean2 += c*c;
    }

    mean /= nc;
    mean2 /= nc;

    double var = mean2 - mean*mean;

    EXPECT_NEAR(1./12, var, 1E-3);
    EXPECT_NEAR(0.5, mean, 1E-3);

    for (uint n = 0; n < N; ++n)
    {
        uint choice = binarySearchForInterval(n + 0.5, cumsum.memptr(), N);

        EXPECT_EQ(n, choice);

        if (HasFailure())
        {
            break;
        }
    }

    rates(N/2+1) = 0;
    rates(N/2+2) = 0;
    rates(N/2+3) = 0;
    rates(N/2+4) = 0;

    double prev = 0;
    for (uint i = 0; i < N; ++i)
    {
        cumsum(i) = prev + rates(i);
        prev = cumsum(i);
    }

    EXPECT_EQ(N/2, binarySearchAndScan(cumsum.memptr(), N, N/2));

    rates.zeros();

    rates(N/2 - 3) = 1;

    prev = 0;
    for (uint i = 0; i < N; ++i)
    {
        cumsum(i) = prev + rates(i);
        prev = cumsum(i);
    }

    for (uint n = 0; n < nc; ++n)
    {
        EXPECT_EQ(N/2 - 3, chooseFromTotalRate(cumsum, 1));
    }
}

}  // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    rng.initialize(time(nullptr));

    return RUN_ALL_TESTS();
}
