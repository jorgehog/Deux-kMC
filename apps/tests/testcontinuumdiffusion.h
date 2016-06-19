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

        m_L = 10;
        m_W = 10;
        m_alpha = 1.0;
        m_mu = 1;
        m_height = 5 + rng.uniform();
        double maxdt = 0.1;

        m_solver = new SOSSolver(m_L, m_W, m_alpha, m_mu);
        setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, m_L, m_W);
        m_pressureWallEvent = new FixedSurface(*m_solver, m_height);
        m_cdiffusionEvent = new RadialFirstPassage(*m_solver, maxdt, 3, 0.038);

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

    RadialFirstPassage *m_cdiffusionEvent;
};

TEST_F(CDiffTest, scan)
{
    const uint n = 0;

    m_solver->setHeight(1, 1, -1, false);
    m_cdiffusionEvent->insertParticle(1, 1, 0);
    m_solver->setHeight(1, 1, 0, false);

    const double &x0 = m_cdiffusionEvent->particlePositions()(0, n);
    const double &y0 = m_cdiffusionEvent->particlePositions()(1, n);
    const double &z0 = m_cdiffusionEvent->particlePositions()(2, n);

    EXPECT_EQ(1, x0);
    EXPECT_EQ(1, y0);
    EXPECT_EQ(0, z0);

    EXPECT_TRUE(solver().isBlockedPosition(x0, y0, z0));

    m_cdiffusionEvent->scan(n, 2, 0.01);

    EXPECT_FALSE(solver().isBlockedPosition(x0, y0, z0)) << m_cdiffusionEvent->particlePositions();

    EXPECT_NEAR(1, x0, 1E-2);
    EXPECT_NEAR(1, y0, 1E-2);
    EXPECT_NEAR(1, z0, 1E-2);

    if (HasFailure())
    {
        cout << x0 << " " << y0 << " " << z0 << endl;
        cout << solver().isBlockedPosition(x0, y0, z0) << endl;
        cout << solver().height(x0, y0) << endl;
        return;
    }

    solver().setHeight(1,1,1,false);

    EXPECT_TRUE(solver().isBlockedPosition(x0, y0, z0));

    m_cdiffusionEvent->scan(n, 0, -0.01, 10000);

    EXPECT_FALSE(solver().isBlockedPosition(x0, y0, z0)) << m_cdiffusionEvent->particlePositions();

    EXPECT_NEAR(0.5, x0, 1E-2);
    EXPECT_NEAR(1, y0, 1E-2);
    EXPECT_NEAR(1, z0, 1E-2);



}


TEST_F(CDiffTest, fullscan)
{
    const uint n = 0;

    m_cdiffusionEvent->insertParticle(1, 1, 1.8);
    m_solver->setHeight(1, 1, 1, false);

    uint dim;
    double delta;
    m_cdiffusionEvent->scanForDisplacement(n, dim, delta);

    EXPECT_EQ(2, dim);
    EXPECT_NEAR(0.2, delta, 1E-2);

    m_cdiffusionEvent->clearDiffusingParticles();

    m_solver->setHeight(1, 1, 0, false);
    m_cdiffusionEvent->insertParticle(1, 0.9, 1);
    m_solver->setHeight(1, 1, 1, false);

    m_cdiffusionEvent->scanForDisplacement(n, dim, delta, 0.01);

    EXPECT_EQ(1, dim);
    EXPECT_NEAR(-0.4, delta, 1E-2);

}


TEST_F(CDiffTest, dissolutionDeposition)
{
    EXPECT_EQ(0, m_cdiffusionEvent->nOfflatticeParticles());

    m_solver->registerHeightChange(1, 1, -1);
    m_cdiffusionEvent->calculateLocalRatesAndUpdateDepositionRates();

    EXPECT_EQ(1, m_cdiffusionEvent->nOfflatticeParticles());

    EXPECT_NEAR(m_cdiffusionEvent->c(), m_cdiffusionEvent->localRates(1, 1, 0), 1E-3);

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
            m_cdiffusionEvent->calculateLocalRatesAndUpdateDepositionRates();

            EXPECT_NEAR(m_cdiffusionEvent->c(), m_cdiffusionEvent->localRates(x, y, m_cdiffusionEvent->nOfflatticeParticles()-1), 1E-3);

            if (HasFailure())
            {
                return;
            }
        }
    }

    m_cdiffusionEvent->releaseLockedParticles();

    EXPECT_EQ(solver().area(), m_cdiffusionEvent->nOfflatticeParticles());

    double T = 100.0;
    m_cdiffusionEvent->diffuseFull(T);
    m_cdiffusionEvent->calculateLocalRatesAndUpdateDepositionRates();

    for (uint x = 0; x < m_L; ++x)
    {
        for (uint y = 0; y < m_W; ++y)
        {
            for (uint n = 0; n < m_cdiffusionEvent->nOfflatticeParticles(); ++n)
            {
                EXPECT_NEAR(m_cdiffusionEvent->_localRateOverD(x, y, n),
                            m_cdiffusionEvent->localRates(x, y, n), 1E-3);

                if (HasFailure())
                {
                    return;
                }
            }
        }
    }

}

TEST_F(CDiffTest, volume)
{
//    solver().setHeights(randi<imat>(m_L, m_W, distr_param(-4, 4)), false);
    solver().confiningSurfaceEvent().setHeight(5);

    const uint N = 10000000;

    double inside = 0;
    const int zMin = solver().heights().min();
    const double zSpan = solver().confiningSurfaceEvent().height() - zMin;

    for (uint n = 0; n < N; ++n)
    {
        const double x = rng.uniform()*m_L - 0.5;
        const double y = rng.uniform()*m_W - 0.5;
        const double z = zMin + rng.uniform()*zSpan;

        if (!solver().isBlockedPosition(x, y, z))
        {
            inside++;
        }
    }

    const double frac = inside/N;

    const double fullBoxVolume = m_L*m_W*zSpan;

    const double soluteVolume = solver().freeVolume();

    const double estimatedVolume = frac*fullBoxVolume;

    const double relErr = fabs(estimatedVolume - soluteVolume)/soluteVolume;

    EXPECT_NEAR(soluteVolume, estimatedVolume, soluteVolume/100.) << soluteVolume << " " << frac*fullBoxVolume;
    EXPECT_NEAR(0, relErr, 1E-3);

    cout << "err: " << relErr << endl;
}























