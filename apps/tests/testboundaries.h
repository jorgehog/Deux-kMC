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
        x0 = getBoundaryFromID(m_solver, x0ID, 0, L, W, Boundary::orientations::FIRST);

        for (uint x1ID = 0; x1ID < isBlocked.size(); ++x1ID)
        {
            x1 = getBoundaryFromID(m_solver, x1ID, 0, L, W, Boundary::orientations::LAST);

            for (uint y0ID = 0; y0ID < isBlocked.size(); ++y0ID)
            {
                y0 = getBoundaryFromID(m_solver, y0ID, 1, W, L, Boundary::orientations::FIRST);

                for (uint y1ID = 0; y1ID < isBlocked.size(); ++y1ID)
                {
                    y1 = getBoundaryFromID(m_solver, y1ID, 1, W, L, Boundary::orientations::LAST);

                    EXPECT_EQ(x0->isBlocked(-1, 0, 0), isBlocked.at(x0ID));
                    EXPECT_EQ(x1->isBlocked(L, 0, 0), isBlocked.at(x1ID));
                    EXPECT_EQ(y0->isBlocked(-1, 0, 0), isBlocked.at(y0ID));
                    EXPECT_EQ(y1->isBlocked(W, 0, 0), isBlocked.at(y1ID));

                    //towards center is never blocked
                    EXPECT_FALSE(x0->isBlocked(1, 1, 1));
                    EXPECT_FALSE(x1->isBlocked(1, 1, 1));
                    EXPECT_FALSE(y0->isBlocked(1, 1, 1));
                    EXPECT_FALSE(y1->isBlocked(1, 1, 1));

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

//--gtest_filter=SOSkMCTest.boundaries_concentration
TEST_F(SOSkMCTest, boundaries_concentration)
{
    const uint L = 10;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 10.23123;
    const int iheight = (int)height;

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

    ConcentrationBoundaryReaction r(0, 0, *m_solver);

    double origArea = W*(height-1);
    EXPECT_EQ(origArea, r.freeBoundaryArea());

    uint yloc = 1;
    int zloc = 5;
    SOSDiffusionReaction *dr = diffusionEvent->addDiffusionReactant(0, yloc, zloc);

    EXPECT_EQ(origArea-1, r.freeBoundaryArea());

    uint n = 0;
    int shift = 0;
    uint yFree;
    int zFree;
    for (uint y = 0; y < W; ++y)
    {
        for (int z = m_solver->height(0, y) + 1; z < iheight; ++z)
        {
            if (y == yloc && z == zloc)
            {
                shift = 1;
            }

            r.getFreeBoundarSite(n, yFree, zFree);

            if (shift == 0)
            {
                EXPECT_EQ(y, yFree);
                EXPECT_EQ(z, zFree);
            }

            else
            {

                if (z == iheight - 1)
                {
                    EXPECT_EQ(y+1, yFree);
                    EXPECT_EQ(m_solver->height(0, y+1) + 1, zFree);
                }

                else
                {
                    EXPECT_EQ(y, yFree);
                    EXPECT_EQ(z+1, zFree);
                }
            }


            n++;

            if (n == r.freeBoundarySites())
            {
                break;
            }

            if (HasFailure())
            {
                return;
            }
        }
    }


    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            for (int z = m_solver->height(x, y)+1; z < iheight - 1; ++z)
            {
                if (x == 0)
                {
                    EXPECT_TRUE(r.pointIsOnBoundary(x, y));
                }

                else
                {
                    EXPECT_FALSE(r.pointIsOnBoundary(x, y));
                }
            }
        }
    }

    r.calculateRate();
    double rate = r.rate();
    for (uint y = 0; y < W; ++y)
    {
        for (int z = m_solver->height(0, y)+2; z < iheight-1; ++z)
        {
            diffusionEvent->executeDiffusionReaction(dr, 0, y, z);
        }

        EXPECT_NEAR(rate, r.rateExpression(), 1E-3);
    }

    diffusionEvent->clearDiffusionReactions();

    n = 0;
    uint max = W*(iheight-1);
    for (uint y = 0; y < W; ++y)
    {
        for (int z = m_solver->height(0, y)+2; z < iheight; ++z)
        {
            diffusionEvent->addDiffusionReactant(0, y, z);
            n++;
        }

        EXPECT_EQ(max-n, r.freeBoundarySites());
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














