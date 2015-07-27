#pragma once

#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

#include "../../src/soskmc/sosdiffusionreaction.h"

TEST_F(SOSkMCTest, diffusion)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double dt = 1.;
    const double height = 20 + rng.uniform();
    const uint spacing = 10;


    Boundary* xBoundary = new Periodic(L);
    Boundary* yBoundary = new Periodic(W);
    m_solver = new SOSSolver(L, W, alpha, mu, xBoundary, yBoundary);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    m_diffusionEvent = new OfflatticeMonteCarloBoundary(*m_solver, dt, spacing);

    SetUp_yo();

    //    OfflatticeMonteCarloBoundary *diffusionEvent = static_cast<OfflatticeMonteCarloBoundary*>(m_diffusionEvent);

    primeSolver(0);

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_EQ(m_solver->nNeighbors(x, y) + m_solver->numberOfSurroundingSolutionSites(x, y), 6);
            m_solver->setHeight(x, y, 0);
        }
    }

    const uint cx = L/2;
    const uint cy = W/2;

    int dx, dy, dz;


    //COMPLETELY FLAT

    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);

    EXPECT_EQ(0, dx); EXPECT_EQ(0, dy); EXPECT_EQ(2, dz);

    //NO NEIGHBORS

    m_solver->registerHeightChange(cx, cy, 1);

    EXPECT_EQ(5, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    //0 = above => dr = (0, 0, 2)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(2, dz);

    //1 = left => dr = (-1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ(-1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //2 = right => dr = (1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //3 = bottom => dr = (0, -1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 3);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(1, dz);

    //4 = top => dr = (0, 1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 4);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 1, dy); EXPECT_EQ(1, dz);


    //ONE NEIGHBOR AT LEFT
    //increase left site to 1
    m_solver->registerHeightChange(cx-1, cy, 1);

    EXPECT_EQ(4, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    //0 = above => dr = (0, 0, 2)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(2, dz);

    //1 = left => dr = (-1, 0, 1) BLOCKED

    //1 = right => dr = (1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ( 1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //2 = bottom => dr = (0, -1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(1, dz);

    //3 = top => dr = (0, 1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 3);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 1, dy); EXPECT_EQ(1, dz);



    //ONE NEIGHBOR AT LEFT AND ONE AT TOP
    //increase bottom site to 2 and top site to 1
    m_solver->registerHeightChange(cx-1, cy, 1);
    m_solver->registerHeightChange(cx, cy+1, 1);

    EXPECT_EQ(3, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    //0 = above => dr = (0, 0, 2)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(2, dz);

    //1 = left => dr = (-1, 0, 1) BLOCKED

    //1 = right => dr = (1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ( 1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //2 = bottom => dr = (0, -1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(1, dz);

    //3 = top => dr = (0, 1, 1) //BLOCKED



    //ONE NEIGHBOR AT RIGHT AND ONE AT TOP
    //increase site to 2 and right site to 2
    m_solver->registerHeightChange(cx, cy, 1);      //increase site to 1
    //2 2 0 0 1
    EXPECT_EQ(4, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    m_solver->registerHeightChange(cx+1, cy, 1);
    m_solver->registerHeightChange(cx+1, cy, 1);    //increase right site to 2
    //2 2 2 0 1
    EXPECT_EQ(3, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    m_solver->registerHeightChange(cx, cy+1, +1);   //increase top site to 2
    //2 2 2 0 2
    EXPECT_EQ(2, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    m_solver->registerHeightChange(cx-1, cy, -1);   //decrease left site to 1
    //1 2 2 0 2
    EXPECT_EQ(3, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    //0 = above => dr = (0, 0, 2)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(2, dz);

    //1 = left => dr = (-1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ(-1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //2 = right => dr = (1, 0, 1) //BLOCKED

    //2 = bottom => dr = (0, -1, 1) //BLOCKED
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(1, dz);

    //2 = top => dr = (0, 1, 1) //BLOCKED


    /*  0        2        0
        1        2        2
        0        0        0  */
}

TEST_F(SOSkMCTest, surfaceSites)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double dt = 1.;
    const double height = 20 + rng.uniform();
    const uint spacing = 10;


    Boundary* xBoundary = new Periodic(L);
    Boundary* yBoundary = new Periodic(W);
    m_solver = new SOSSolver(L, W, alpha, mu, xBoundary, yBoundary);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    OfflatticeMonteCarloBoundary *diffusionEvent = new OfflatticeMonteCarloBoundary(*m_solver, dt, spacing);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();
    rng.initialize(1000);

    //    OfflatticeMonteCarloBoundary *diffusionEvent = static_cast<OfflatticeMonteCarloBoundary*>(m_diffusionEvent);

    primeSolver(0);


    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_EQ(m_solver->nNeighbors(x, y) + m_solver->numberOfSurroundingSolutionSites(x, y), 6);
            m_solver->setHeight(x, y, 0);
        }
    }

    diffusionEvent->clearDiffusionReactions();

    EXPECT_EQ(0, diffusionEvent->numberOfDiffusionReactions());

    uint c = 0;
    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            SOSDiffusionReaction *r = diffusionEvent->addDiffusionReactant(x, y, 2);
            EXPECT_EQ(x, r->x()); EXPECT_EQ(y, r->y()); EXPECT_EQ(2, r->z());

            c++;
            EXPECT_EQ(c, diffusionEvent->numberOfDiffusionReactions()) << "wrong number of diffusionreactions added: " << x << " " << y;
        }
    }

    EXPECT_EQ(9, diffusionEvent->numberOfDiffusionReactions());


    /*  0        0        0
        0        0        0
        0        0        0  */

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_TRUE(m_solver->isSurfaceSite(x, y, 1));
        }
    }

    SOSDiffusionReaction *r = diffusionEvent->diffusionReaction(1, 1, 2);
    EXPECT_EQ(6, r->numberOfFreePaths());

    //Deposit 4 paticles (neighbors of center site)
    diffusionEvent->diffusionReaction(0, 1, 2)->executeReaction(0, 0, -1);
    diffusionEvent->diffusionReaction(2, 1, 2)->executeReaction(0, 0, -1);
    diffusionEvent->diffusionReaction(1, 0, 2)->executeReaction(0, 0, -1);
    diffusionEvent->diffusionReaction(1, 2, 2)->executeReaction(0, 0, -1);
    EXPECT_EQ(6, r->numberOfFreePaths());

    EXPECT_EQ(1, m_solver->height(0, 1));
    EXPECT_EQ(1, m_solver->height(2, 1));
    EXPECT_EQ(1, m_solver->height(1, 0));
    EXPECT_EQ(1, m_solver->height(1, 2));

    //Deposit 4 particles (not neighbors of center, to neighbors of center)
    diffusionEvent->diffusionReaction(0, 0, 2)->executeReaction( 1, 0, 0);
    EXPECT_EQ(5, r->numberOfFreePaths());

    diffusionEvent->diffusionReaction(2, 2, 2)->executeReaction(-1, 0, 0);
    EXPECT_EQ(4, r->numberOfFreePaths());

    diffusionEvent->diffusionReaction(0, 2, 2)->executeReaction( 0,-1, 0);
    EXPECT_EQ(3, r->numberOfFreePaths());

    diffusionEvent->diffusionReaction(2, 0, 2)->executeReaction( 0, 1, 0);
    EXPECT_EQ(2, r->numberOfFreePaths());

    EXPECT_EQ(2, m_solver->height(0, 1));
    EXPECT_EQ(2, m_solver->height(2, 1));
    EXPECT_EQ(2, m_solver->height(1, 0));
    EXPECT_EQ(2, m_solver->height(1, 2));

    EXPECT_EQ(1, diffusionEvent->numberOfDiffusionReactions());

    r->setZ((uint)height);
    EXPECT_EQ(5, r->numberOfFreePaths());

    r->setZ(2);
    EXPECT_EQ(2, r->numberOfFreePaths());

    uint upCount = 0;
    uint downCount = 0;
    const uint N = 1000000;

    int dx, dy, dz;

    for (uint i = 0; i < N; ++i)
    {
        r->getRandomDiffusionPath(dx, dy, dz);

        EXPECT_EQ(0, dx);
        EXPECT_EQ(0, dy);

        if (dx != 0 || dy != 0)
        {
            break;
        }

        if (dz == 1)
        {
            upCount += 1;
        }

        else if (dz == -1)
        {
            downCount += 1;
        }

        else
        {
            EXPECT_EQ(0, 1);
        }
    }

    double upProb = upCount/double(N);
    double downProb = downCount/double(N);

    EXPECT_NEAR(0.5, upProb, 1E-3);
    EXPECT_NEAR(0.5, downProb, 1E-3);

    r->executeReaction(0, 0, -1);

    EXPECT_EQ(0, diffusionEvent->numberOfDiffusionReactions());

}
