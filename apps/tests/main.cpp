#include <kMC.h>
#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

namespace {

// The fixture for testing class kMCTest.
class SOSkMCTest : public ::testing::Test
{
protected:
    // You can remove any or all of the following functions if its body
    // is empty.

    SOSkMCTest()
    {
        // You can do set-up work for each test here.
    }

    virtual ~SOSkMCTest()
    {
        // You can do clean-up work that doesn't throw exceptions here.
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp()
    {
        // Code here will be called immediately after the constructor (right
        // before each test).

        uint L = 10;
        uint W = 10;
        double alpha = 1.0;
        double mu = 0;
        double E0 = 0.01;
        double sigma0 = 1;
        double r0 = 4;
        double D = 10.;
        double dt = 1.;

        Boundary* xBoundary = new Periodic(L);
        Boundary* yBoundary = new Periodic(W);
        m_solver = new SolidOnSolidSolver(L, W, xBoundary, yBoundary, alpha, mu);
        m_pressureWallEvent = new PressureWall(*m_solver, E0, sigma0, r0);
        m_diffusionEvent = new CavityDiffusion(*m_solver, D, dt);

        m_averageHeight = new AverageHeight(*m_solver);

        m_pressureWallEvent->setDependency(m_averageHeight);

        m_lattice = new Lattice();

        m_lattice->addEvent(m_solver);
        m_lattice->addEvent(m_averageHeight);
        m_lattice->addEvent(m_pressureWallEvent);
        m_lattice->addEvent(m_diffusionEvent);

        m_lattice->enableOutput(false);
        m_lattice->enableEventValueStorage(false, false);

    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test (right
        // before the destructor).

        delete m_solver;
        delete m_pressureWallEvent;
        delete m_averageHeight;
        delete m_diffusionEvent;

        delete m_lattice;
    }

    // Objects declared here can be used by all tests in the test case for Foo.

    void primeSolver(const uint N)
    {
        BasicEvent<uint> primer("kmc primer", [&] (BasicEvent<uint> *event)
        {
            if (event->cycle() == N)
            {
                event->stopLoop();
            }
        });

        m_lattice->addEvent(primer);

        m_lattice->eventLoop(2*N);
    }

    SolidOnSolidSolver *m_solver;
    AverageHeight *m_averageHeight;
    PressureWall *m_pressureWallEvent;
    CavityDiffusion *m_diffusionEvent;

    Lattice *m_lattice;
};

TEST_F(SOSkMCTest, test)
{
    uint prime = 10;
    primeSolver(prime);

    uint repCount = 100;
    for (uint rep = 0; rep < repCount; ++rep)
    {

        uint x = rng.uniform()*m_solver->length();
        uint y = rng.uniform()*m_solver->width();

        const int originalHeight = m_solver->height(x, y);
        const uint originalNParticles = m_diffusionEvent->nParticles();

        const mat originalParticlePositions = m_diffusionEvent->particlePositions();

        m_solver->registerHeightChange(x, y, -1);

        const int newHeight = m_solver->height(x, y);
        const uint newNParticles = m_diffusionEvent->nParticles();
        const mat newParticlePositions = m_diffusionEvent->particlePositions();

        EXPECT_EQ(newHeight, originalHeight - 1);
        EXPECT_EQ(newNParticles, originalNParticles + 1);

        for (uint n = 0; n < originalNParticles; ++n)
        {
            if (originalParticlePositions(2, n) >= m_pressureWallEvent->height() - 0.5)
            {
                continue;
            }

            EXPECT_EQ(newParticlePositions(0, n), originalParticlePositions(0, n)) << "X fail at n=" << n << "/" << newNParticles << " rep = " << rep;
            EXPECT_EQ(newParticlePositions(1, n), originalParticlePositions(1, n)) << "Y fail at n=" << n << "/" << newNParticles << " rep = " << rep;
            EXPECT_EQ(newParticlePositions(2, n), originalParticlePositions(2, n)) << "Z fail at n=" << n << "/" << newNParticles << " rep = " << rep;
        }

        double xNew = m_diffusionEvent->particlePositions()(0, originalNParticles);
        double yNew = m_diffusionEvent->particlePositions()(1, originalNParticles);
        double zNew = m_diffusionEvent->particlePositions()(2, originalNParticles);

        double dr = pow(xNew - (double)x, 2) + pow(yNew - (double)y, 2) + pow(zNew - (double)originalHeight, 2);

        EXPECT_NEAR(dr, 1, 1E-10);

        m_solver->registerHeightChange(x, y, +1);

        const int newHeight2 = m_solver->height(x, y);
        const uint newNParticles2 = m_diffusionEvent->nParticles();
        const mat newParticlePositions2 = m_diffusionEvent->particlePositions();

        EXPECT_EQ(newHeight2, originalHeight);
        EXPECT_EQ(newNParticles2, originalNParticles);

        uint match = 0;
        for (uint n1 = 0; n1 < originalNParticles; ++n1)
        {
            double x0 = originalParticlePositions(0, n1);
            double y0 = originalParticlePositions(1, n1);
            double z0 = originalParticlePositions(2, n1);

            for (uint n2 = 0; n2 < originalNParticles; ++n2)
            {
                double x1 = newParticlePositions2(0, n2);
                double y1 = newParticlePositions2(1, n2);
                double z1 = newParticlePositions2(2, n2);

                if ((x1 == x0) && (y1 == y0) && (z0 == z1))
                {
                    match++;
                    break;
                }
            }
        }

        EXPECT_GE(match, originalNParticles - 1);

        uint nDiff = 100;
        for (uint i = 0; i < nDiff; ++i)
        {
            m_diffusionEvent->diffuse(m_diffusionEvent->dt());
        }

        m_diffusionEvent->_dump(prime + rep);
    }

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
