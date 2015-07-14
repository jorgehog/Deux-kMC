#include <kMC.h>
#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

#include "testdiffusion.h"

namespace {

//TEST_F(SOSkMCTest, diffusionTest)
//{
//    uint prime = 10;
//    primeSolver(prime);

//    uint nDiff = 1000;

//    uint nDump = 0;
//    for (uint n = 0; n < nDiff; ++n)
//    {
//        m_diffusionEvent->diffuse(m_diffusionEvent->dt());

//        cout.flush();

//        if (n % 1 == 0)
//        {
//            m_diffusionEvent->dump(nDump);
//            cout << "\rProgress: " << n/double(nDiff)*100 << " % - acceptance: " << m_diffusionEvent->acceptanceRatio()*100 << " %";
//            nDump++;
//        }
//    }
//    cout << endl;

//}

//TEST_F(SOSkMCTest, test)
//{
//    uint prime = 10;
//    primeSolver(prime);

//    uint repCount = 100;
//    for (uint rep = 0; rep < repCount; ++rep)
//    {

//        uint x = rng.uniform()*m_solver->length();
//        uint y = rng.uniform()*m_solver->width();

//        const int originalHeight = m_solver->height(x, y);
//        const uint originalNParticles = m_diffusionEvent->nParticles();

//        const mat originalParticlePositions = m_diffusionEvent->particlePositions();

//        m_solver->registerHeightChange(x, y, -1);

//        const int newHeight = m_solver->height(x, y);
//        const uint newNParticles = m_diffusionEvent->nParticles();
//        const mat newParticlePositions = m_diffusionEvent->particlePositions();

//        EXPECT_EQ(newHeight, originalHeight - 1);
//        EXPECT_EQ(newNParticles, originalNParticles + 1);


//        for (uint n = 0; n < originalNParticles; ++n)
//        {
//            if (originalParticlePositions(2, n) >= m_pressureWallEvent->height() - 0.5)
//            {
//                continue;
//            }

//            EXPECT_EQ(newParticlePositions(0, n), originalParticlePositions(0, n)) << "X fail at n=" << n << "/" << newNParticles << " rep = " << rep;
//            EXPECT_EQ(newParticlePositions(1, n), originalParticlePositions(1, n)) << "Y fail at n=" << n << "/" << newNParticles << " rep = " << rep;
//            EXPECT_EQ(newParticlePositions(2, n), originalParticlePositions(2, n)) << "Z fail at n=" << n << "/" << newNParticles << " rep = " << rep;
//        }

//        double xNew = m_diffusionEvent->particlePositions()(0, originalNParticles);
//        double yNew = m_diffusionEvent->particlePositions()(1, originalNParticles);
//        double zNew = m_diffusionEvent->particlePositions()(2, originalNParticles);

//        double dr = pow(xNew - (double)x, 2) + pow(yNew - (double)y, 2) + pow(zNew - (double)originalHeight, 2);

//        EXPECT_NEAR(dr, 1, 1E-10);

//        m_solver->registerHeightChange(x, y, +1);

//        const int newHeight2 = m_solver->height(x, y);
//        const uint newNParticles2 = m_diffusionEvent->nParticles();
//        const mat newParticlePositions2 = m_diffusionEvent->particlePositions();

//        EXPECT_EQ(newHeight2, originalHeight);
//        EXPECT_EQ(newNParticles2, originalNParticles);

//        uint match = 0;
//        int sub = 0;
//        for (uint n1 = 0; n1 < originalNParticles; ++n1)
//        {
//            double x0 = newParticlePositions(0, n1);
//            double y0 = newParticlePositions(1, n1);
//            double z0 = newParticlePositions(2, n1);

//            uint XX = newParticlePositions(0, n1);
//            uint YY = newParticlePositions(1, n1);

//            if (XX == x && YY == y && (newParticlePositions(2, n1) < newHeight2 + 0.5))
//            {
//                sub--;
//            }

//            for (uint n2 = 0; n2 < originalNParticles; ++n2)
//            {
//                double x1 = newParticlePositions2(0, n2);
//                double y1 = newParticlePositions2(1, n2);
//                double z1 = newParticlePositions2(2, n2);

//                if ((x1 == x0) && (y1 == y0) && (z0 == z1))
//                {
//                    match++;
//                    break;
//                }
//            }
//        }

//        EXPECT_GE(match, originalNParticles + sub - 1) << match <<" "<< originalNParticles << " " << sub << endl;

//        uint nDiff = 100;
//        for (uint i = 0; i < nDiff; ++i)
//        {
//            m_diffusionEvent->diffuse(m_diffusionEvent->dt());
//        }

//    }

//}

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
