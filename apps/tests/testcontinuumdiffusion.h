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

        m_L = 3;
        m_W = 3;
        m_alpha = 1.0;
        m_mu = 1;
        m_height = 20 + rng.uniform();
        double maxdt = 0.01;
        auto rf = [] (const FirstPassageContinuum *_this, const uint x, const uint y, const uint n)
        {
            return 1.0/_this->nOfflatticeParticles();

            const int z = _this->solver().height(x, y);

            const double dx2 = pow(x - _this->particlePositions(0, n), 2);
            const double dy2 = pow(y - _this->particlePositions(1, n), 2);
            const double dz2 = pow(z - _this->particlePositions(2, n), 2);

            const double dr2 = dx2 + dy2 + dz2;

            return 1./dr2;
        };

        m_solver = new SOSSolver(m_L, m_W, m_alpha, m_mu, getBoundariesFromIDs({0, 0, 0, 0}, m_L, m_W));
        m_pressureWallEvent = new FixedSurface(*m_solver, m_height);
        m_cdiffusionEvent = new FirstPassageContinuum(*m_solver, maxdt, rf);

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

    const double x0Scan = x0;
    const double y0Scan = y0;
    const double z0Scan = z0 + 1;

    m_cdiffusionEvent->scan(n, 2, 0.01);

    EXPECT_FALSE(solver().isBlockedPosition(x0, y0, z0)) << m_cdiffusionEvent->particlePositions();

    EXPECT_NEAR(x0, x0Scan, 1E-3);
    EXPECT_NEAR(y0, y0Scan, 1E-3);
    EXPECT_NEAR(z0, z0Scan, 1E-3);

    solver().setHeight(1,1,1,false);

    EXPECT_TRUE(solver().isBlockedPosition(x0, y0, z0));

    m_cdiffusionEvent->scan(n, 0, -0.01);

    EXPECT_FALSE(solver().isBlockedPosition(x0, y0, z0)) << m_cdiffusionEvent->particlePositions();

    EXPECT_NEAR(0.5, x0, 1E-3);
    EXPECT_NEAR(y0, y0Scan, 1E-3);
    EXPECT_NEAR(z0, z0Scan, 1E-3);



}
































