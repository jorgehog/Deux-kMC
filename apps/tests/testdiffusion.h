#pragma once

#include <kMC.h>
#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

TEST_F(SOSkMCTest, diffusion)
{
    const uint L = 10;
    const uint W = 10;
    const double alpha = 1.0;
    const double mu = 0;
    const double dt = 1.;
    const double height = 20 + rng.uniform();
    const uint spacing = 3;


    Boundary* xBoundary = new Periodic(L);
    Boundary* yBoundary = new Periodic(W);
    m_solver = new SOSSolver(L, W, alpha, mu, xBoundary, yBoundary);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    m_diffusionEvent = new OfflatticeMonteCarloBoundary(*m_solver, dt, spacing);

    SetUp_yo();

    OfflatticeMonteCarloBoundary *diffusionEvent = static_cast<OfflatticeMonteCarloBoundary*>(m_diffusionEvent);

    primeSolver(0);

    EXPECT_TRUE(diffusionEvent->checkIfEnoughRoom());

    const int h0 = m_solver->height(L/2, W/2);
    const int hmax = m_solver->heights().max();
    const double dh = m_solver->confiningSurfaceEvent().height() - hmax;

    for (int h = h0; h <= (int)(dh - 2*spacing); ++h)
    {
        m_solver->registerHeightChange(L/2, W/2, 1);
        EXPECT_TRUE(diffusionEvent->checkIfEnoughRoom());
    }

    m_solver->registerHeightChange(L/2, W/2, 1);
    EXPECT_FALSE(diffusionEvent->checkIfEnoughRoom());


    m_lattice->reConnect();

}
