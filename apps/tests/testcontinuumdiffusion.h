#pragma once

#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

class CDiffTest : public SOSkMCTest
{
public:
    void SetUp()
    {
        SOSkMCTest::SetUp();

        m_L = 5;
        m_W = 5;
        m_alpha = 1.0;
        m_mu = 1;
        m_height = 5 + rng.uniform();
        double maxdt = 0.01;

        m_solver = new SOSSolver(m_L, m_W, m_alpha, m_mu);
        setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, m_L, m_W);
        m_pressureWallEvent = new FixedSurface(*m_solver, m_height);
        m_cdiffusionEvent = new FirstPassageContinuum(*m_solver, maxdt, 2, 1);

        m_diffusionEvent = m_cdiffusionEvent;
        SetUp_yo();

        for (uint x = 0; x < m_L; ++x)
        {
            for (uint y = 0; y < m_W; ++y)
            {
                m_solver->setHeight(x, y, 0, false);
            }
        }

        primeSolver();

        m_cdiffusionEvent->clearDiffusingParticles();


    }

    uint m_L;
    uint m_W;
    double m_alpha;
    double m_mu;
    double m_height;

    FirstPassageContinuum *m_cdiffusionEvent;
};

TEST_F(CDiffTest, scan)
{
    const uint n = 0;

    m_cdiffusionEvent->insertParticle(1, 1, 0);

    const double &x0 = m_cdiffusionEvent->particlePositions()(0, n);
    const double &y0 = m_cdiffusionEvent->particlePositions()(1, n);
    const double &z0 = m_cdiffusionEvent->particlePositions()(2, n);

    EXPECT_TRUE(solver().isBlockedPosition(x0, y0, z0));

    m_cdiffusionEvent->scan(n, 2, 0.01);

    EXPECT_FALSE(solver().isBlockedPosition(x0, y0, z0)) << m_cdiffusionEvent->particlePositions();

    EXPECT_NEAR(1, x0, 1E-2);
    EXPECT_NEAR(1, y0, 1E-2);
    EXPECT_NEAR(0.5, z0, 1E-2);

    solver().setHeight(1,1,1,false);

    EXPECT_TRUE(solver().isBlockedPosition(x0, y0, z0));

    m_cdiffusionEvent->scan(n, 0, -0.01, 10000);

    EXPECT_FALSE(solver().isBlockedPosition(x0, y0, z0)) << m_cdiffusionEvent->particlePositions();

    EXPECT_NEAR(0.5, x0, 1E-2);
    EXPECT_NEAR(1, y0, 1E-2);
    EXPECT_NEAR(0.5, z0, 1E-2);



}


TEST_F(CDiffTest, fullscan)
{
    const uint n = 0;

    m_solver->setHeight(1, 1, 1, false);

    m_cdiffusionEvent->insertParticle(1, 1, 1.3);

    uint dim;
    double delta;
    m_cdiffusionEvent->scanForDisplacement(n, dim, delta);

    EXPECT_EQ(2, dim);
    EXPECT_NEAR(0.2, delta, 1E-2);

    m_cdiffusionEvent->clearDiffusingParticles();

    m_cdiffusionEvent->insertParticle(1, 0.9, 1);

    m_cdiffusionEvent->scanForDisplacement(n, dim, delta, 0.01);

    EXPECT_EQ(1, dim);
    EXPECT_NEAR(-0.4, delta, 1E-2);

}


TEST_F(CDiffTest, dissolutionDeposition)
{
    EXPECT_EQ(0, m_cdiffusionEvent->nOfflatticeParticles());

    cout << m_solver->closestSquareDistance(1, 1, 0, 1, 1, 1) << endl;

    m_solver->registerHeightChange(1, 1, -1);

    EXPECT_EQ(1, m_cdiffusionEvent->nOfflatticeParticles());

    EXPECT_NEAR(1, m_cdiffusionEvent->localRates(1, 1, 0), 1E-3);

    m_solver->registerHeightChange(1, 1, 1);

    EXPECT_EQ(0, m_cdiffusionEvent->nOfflatticeParticles());
}

TEST_F(CDiffTest, diffusion)
{
    for (uint x = 0; x < m_L; ++x)
    {
        for (uint y = 0; y < m_W; ++y)
        {
            solver().registerHeightChange(x, y, -1);

            EXPECT_NEAR(1, m_cdiffusionEvent->localRates(x, y, m_cdiffusionEvent->nOfflatticeParticles()-1), 1E-3);

            if (HasFailure())
            {
                return;
            }
        }
    }

    EXPECT_EQ(solver().area(), m_cdiffusionEvent->nOfflatticeParticles());

    double T = 100.0;
    m_cdiffusionEvent->diffuseFull(T);

    EXPECT_NEAR(1, m_cdiffusionEvent->acceptanceRatio(), 1E-3);

    const uint N = T/m_cdiffusionEvent->maxdt();

    EXPECT_NEAR(N*solver().area(), m_cdiffusionEvent->trials(), 1E-3);

    for (uint x = 0; x < m_L; ++x)
    {
        for (uint y = 0; y < m_W; ++y)
        {
            for (uint n = 0; n < m_cdiffusionEvent->nOfflatticeParticles(); ++n)
            {
                EXPECT_NEAR(m_cdiffusionEvent->calculateLocalRate(x, y, n),
                            m_cdiffusionEvent->localRates(x, y, n), 1E-3);
            }
        }
    }

}

TEST_F(CDiffTest, volume)
{
    solver().setHeights(randi<imat>(m_L, m_W, distr_param(-4, 4)), false);

    const uint N = 10000000;

    double inside = 0;
    const int zMin = solver().heights().min();
    const double zSpan = solver().confiningSurfaceEvent().height() - zMin;

    const double correction = solver().calculateVolumeCorrection();

    for (uint n = 0; n < N; ++n)
    {
        const double x = rng.uniform()*m_L;
        const double y = rng.uniform()*m_W;
        const double z = zMin + rng.uniform()*zSpan;

        if (!solver().isBlockedPosition(x, y, z))
        {
            inside++;
        }
    }

    const double frac = inside/N;

    const double fullBoxVolume = m_L*m_W*zSpan;

    const double soluteVolume = solver().volume();

    const double estimatedVolume = frac*fullBoxVolume + correction;

    const double relErr = fabs(estimatedVolume - soluteVolume)/soluteVolume;

    EXPECT_NEAR(soluteVolume, estimatedVolume, soluteVolume/100.) << correction << " " << soluteVolume << " " << frac*fullBoxVolume;
    EXPECT_NEAR(0, relErr, 1E-3);

    cout << "err: " << relErr << endl;
}






















