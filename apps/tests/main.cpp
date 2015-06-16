#include <kMC.h>
#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

namespace {

// The fixture for testing class kMCTest.
class kMCTest : public ::testing::Test
{
protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    kMCTest()
    {
        // You can do set-up work for each test here.
    }

    virtual ~kMCTest()
    {
        // You can do clean-up work that doesn't throw exceptions here.
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp()
    {
        // Code here will be called immediately after the constructor (right
        // before each test).
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }

    // Objects declared here can be used by all tests in the test case for Foo.
};

TEST_F(kMCTest, test)
{
    EXPECT_EQ(0, 1-1);
}

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

    EXPECT_NEAR(1./12, var, 1E-3);
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
