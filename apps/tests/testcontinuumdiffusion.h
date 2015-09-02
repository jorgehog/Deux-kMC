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
        m_mu = 0;
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
        FirstPassageContinuum *diffusionEvent = new FirstPassageContinuum(*m_solver, maxdt, rf);

        m_diffusionEvent = diffusionEvent;
        SetUp_yo();

        primeSolver(0);

        diffusionEvent->clearDiffusingParticles();

        for (uint x = 0; x < m_L; ++x)
        {
            for (uint y = 0; y < m_W; ++y)
            {
                m_solver->setHeight(x, y, 0);
                diffusionEvent->clearDiffusingParticles();
            }
        }
    }

private:

    uint m_L;
    uint m_W;
    double m_alpha;
    double m_mu;
    double m_height;
};

TEST_F(CDiffTest, scan)
{

}

