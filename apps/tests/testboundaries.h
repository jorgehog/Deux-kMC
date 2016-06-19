#pragma once

#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

#include "../../src/soskmc/sosdiffusionreaction.h"

#include "../../src/soskmc/concentrationboundaryreaction.h"

//--gtest_filter=SOSkMCTest.boundaries
TEST_F(SOSkMCTest, boundaries_blocked)
{
    //    rng.initialize(1230303);

    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 1, 2}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;
    SetUp_yo();

    primeSolver();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    vector<bool> isBlocked = {false, true, false};

    const Boundary *x0;
    const Boundary *x1;
    const Boundary *y0;
    const Boundary *y1;

    for (uint x0ID = 0; x0ID < isBlocked.size(); ++x0ID)
    {
        x0 = getBoundaryFromID(x0ID, L, Boundary::orientations::FIRST);

        for (uint x1ID = 0; x1ID < isBlocked.size(); ++x1ID)
        {
            x1 = getBoundaryFromID(x1ID, L, Boundary::orientations::LAST);

            for (uint y0ID = 0; y0ID < isBlocked.size(); ++y0ID)
            {
                y0 = getBoundaryFromID(y0ID, W, Boundary::orientations::FIRST);

                for (uint y1ID = 0; y1ID < isBlocked.size(); ++y1ID)
                {
                    y1 = getBoundaryFromID(y1ID, W, Boundary::orientations::LAST);

                    EXPECT_EQ(x0->isBlockedLattice(-1, 0, 0), isBlocked.at(x0ID));
                    EXPECT_EQ(x1->isBlockedLattice(L, 0, 0), isBlocked.at(x1ID));
                    EXPECT_EQ(y0->isBlockedLattice(-1, 0, 0), isBlocked.at(y0ID));
                    EXPECT_EQ(y1->isBlockedLattice(W, 0, 0), isBlocked.at(y1ID));

                    //towards center is never blocked
                    EXPECT_FALSE(x0->isBlockedLattice(1, 1, 1));
                    EXPECT_FALSE(x1->isBlockedLattice(1, 1, 1));
                    EXPECT_FALSE(y0->isBlockedLattice(1, 1, 1));
                    EXPECT_FALSE(y1->isBlockedLattice(1, 1, 1));

                    if (HasFailure())
                    {
                        return;
                    }

                    delete y1;
                }

                delete y0;
            }

            delete x1;
        }

        delete x0;
    }


    EXPECT_EQ(m_solver->calculateNNeighbors(0, 1), 5);
    EXPECT_EQ(m_solver->calculateNNeighbors(2, 1), 5);
    EXPECT_EQ(m_solver->calculateNNeighbors(1, 0), 5);
    EXPECT_EQ(m_solver->calculateNNeighbors(1, 2), 4);

    EXPECT_EQ(m_solver->nNeighbors(0, 1), 5);
    EXPECT_EQ(m_solver->nNeighbors(2, 1), 5);
    EXPECT_EQ(m_solver->nNeighbors(1, 0), 5);
    EXPECT_EQ(m_solver->nNeighbors(1, 2), 4);
}


void testDiffYo(RadialFirstPassage *diff)
{
#ifndef NDEBUG
    return;
#endif

    const uint np = 1000;
    for (uint p = 0; p < np; ++p)
    {
        diff->insertRandomParticle();
    }

    uint nH = 20;
    mat hist(nH, nH, fill::zeros);
    const double delta = const_cast<const RadialFirstPassage*>(diff)->solver().length()/double(nH);

    const uint N = 1000000;

    //    double x, y, z;
    for (uint n = 0; n < N; ++n)
    {
        diff->diffuse(0.01);

        for (uint p = 0; p < np; ++p)
        {
            //            diff->findRandomPosition(x, y, z);
            const double &x = diff->particlePositions(0, p);
            const double &y = diff->particlePositions(1, p);

            const uint nx = uint((x+0.5)/delta);
            const uint ny = uint((y+0.5)/delta);

            hist(nx, ny)++;
        }

        if (n % 100 == 0)
        {
            cout.flush();
            cout << "\r" << n/double(N-1);
        }
    }

    hist /= N;
    hist /= hist.max();

    cout << endl;
    cout << hist << endl;

    for (uint i = 0; i < nH; ++i)
    {
        for (uint j = 0; j < nH; ++j)
        {
            const double hval = hist(i, j);

            EXPECT_NEAR(1, hval, 0.05);
        }
    }
}

TEST_F(SOSkMCTest, boundaries_refl)
{
    rng.initialize(time(nullptr));

    const uint L = 30;
    const uint W = 30;
    const double alpha = 3.0;
    const double gamma = 0;
    const double height = 10.3;

    m_solver = new SOSSolver(L, W, alpha, gamma, true);
    setBoundariesFromIDs(m_solver, {3, 3, 3, 3}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    RadialFirstPassage *diff = new RadialFirstPassage(*m_solver, 0.01, 3, 0.1);
    m_diffusionEvent = diff;

    initializeSurface(*m_solver, "flat");

    SetUp_yo();
    primeSolver();
    diff->clearDiffusingParticles();

    testDiffYo(diff);
}


TEST_F(SOSkMCTest, boundaries_periodic)
{
    rng.initialize(time(nullptr));

    const uint L = 30;
    const uint W = 30;
    const double alpha = 3.0;
    const double gamma = 0;
    const double height = 10.3;

    m_solver = new SOSSolver(L, W, alpha, gamma, true);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    RadialFirstPassage *diff = new RadialFirstPassage(*m_solver, 0.01, 3, 0.1);
    m_diffusionEvent = diff;

    initializeSurface(*m_solver, "flat");

    SetUp_yo();
    primeSolver();
    diff->clearDiffusingParticles();

    testDiffYo(diff);
}


//--gtest_filter=SOSkMCTest.boundaries_concentration
TEST_F(SOSkMCTest, boundaries_concentration)
{
    const uint L = 10;
    const uint W = 10;
    const double alpha = 1.0;
    const double gamma = 0;
    const double height = 10.3;
    const int iheight = (int)height;

    m_solver = new SOSSolver(L, W, alpha, gamma, true);
    setBoundariesFromIDs(m_solver, {0, 0, 2, 2}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    RadialFirstPassage *diff = new RadialFirstPassage(*m_solver, 0.01, 3, 0.1);
    m_diffusionEvent = diff;

    initializeSurface(*m_solver, "flat");

    SetUp_yo();
    primeSolver();

    ConcentrationBoundaryReaction r(0, 0, *m_solver, 0);

    EXPECT_EQ(W*(iheight-1), r.nBoundarySites());

    uint xi;
    int h;
    for (uint n = 0; n < r.nBoundarySites(); ++n)
    {
        r.getBoundaryPosition(xi, h, n, iheight);

        EXPECT_EQ(n/(iheight - 1), xi);
    }

    ConcentrationBoundaryReaction r2(0, 1, *m_solver, 0);

    EXPECT_EQ(W*(iheight-1), r2.nBoundarySites());

    for (uint n = 0; n < r2.nBoundarySites(); ++n)
    {
        r2.getBoundaryPosition(xi, h, n, iheight);

        EXPECT_EQ(n/(iheight - 1), xi);
    }

    ConcentrationBoundaryReaction r3(1, 0, *m_solver, 0);

    EXPECT_EQ(L*(iheight-1), r3.nBoundarySites());

    for (uint n = 0; n < r3.nBoundarySites(); ++n)
    {
        r3.getBoundaryPosition(xi, h, n, iheight);

        EXPECT_EQ(n/(iheight - 1), xi);
    }

    ConcentrationBoundaryReaction r4(1, 1, *m_solver, 0);

    EXPECT_EQ(L*(iheight-1), r4.nBoundarySites());

    for (uint n = 0; n < r4.nBoundarySites(); ++n)
    {
        r4.getBoundaryPosition(xi, h, n, iheight);

        EXPECT_EQ(n/(iheight - 1), xi);
    }


}


//--gtest_filter=SOSkMCTest.diffusion
TEST_F(SOSkMCTest, boundaries_reflect)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 10;

    const uint bType = 3;

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {bType, bType, bType, bType}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;
    SetUp_yo();

    primeSolver();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    EXPECT_EQ(2, solver().rightSite(2, 1, 0));
    EXPECT_EQ(0, solver().leftSite(0, 1, 0));
    EXPECT_EQ(2, solver().topSite(1, 2, 0));
    EXPECT_EQ(0, solver().bottomSite(1, 0, 0));

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_EQ(5, m_solver->calculateNNeighbors(x, y)) << x << " " << y;
        }
    }

    SOSDiffusionReaction *r = diffusionEvent->addDiffusionReactant(2, 2, 8);

    EXPECT_EQ(4, r->calculateNumberOfFreePaths());

    solver().registerHeightChange(0, 0, 1);
    solver().registerHeightChange(0, 0, 1);
    EXPECT_EQ(3, solver().numberOfSurroundingSolutionSites(0, 0));

}


TEST_F(SOSkMCTest, boundaries_constant_height)
{
    return;

    //THIS TEST IS DEPRECATED


    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 0;

    const uint bType = 4;

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {bType, bType, bType, bType}, L, W, height);
    m_pressureWallEvent = new NoConfinement(*m_solver);
    m_diffusionEvent = new ConstantConcentration(*m_solver);
    SetUp_yo();

    primeSolver();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
        }
    }

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_EQ(5, solver().nNeighbors(x, y));
        }
    }

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->registerHeightChange(x, y, 1);
        }
    }

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_EQ(3 + x%(L-1) + y%(W-1), solver().nNeighbors(x, y));
        }
    }

}














