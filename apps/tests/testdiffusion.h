#pragma once

#include <SOSkMC.h>
#include <gtest/gtest.h>

#include <sys/time.h>

#include "kmctestfixture.h"

#include "../../src/soskmc/sosdiffusionreaction.h"

#include "../../src/soskmc/dissolutiondeposition.h"

class Dummy : public Reaction
{

    // Reaction interface
public:
    bool isAllowed() const {return true;}
    void executeAndUpdate() {}
    double rateExpression() {return 1.0;}
};

TEST_F(SOSkMCTest, diff_reac_map)
{

}

TEST_F(SOSkMCTest, affected_set)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent; SetUp_yo();

    Dummy r;

    m_solver->registerAffectedReaction(&r);

    EXPECT_EQ(1, solver().affectedReactions().size());

    m_solver->registerAffectedReaction(&r);

    EXPECT_EQ(1, solver().affectedReactions().size());

    Dummy r2;

    m_solver->registerAffectedReaction(&r2);

    EXPECT_EQ(2, solver().affectedReactions().size());

    m_solver->registerAffectedReaction(&r2);

    EXPECT_EQ(2, solver().affectedReactions().size());



}

TEST_F(SOSkMCTest, diffusion)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;
    SetUp_yo();

    //    OfflatticeMonteCarloBoundary *diffusionEvent = static_cast<OfflatticeMonteCarloBoundary*>(m_diffusionEvent);

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            uint nn = m_solver->nNeighbors(x, y);
            uint ns = m_solver->numberOfSurroundingSolutionSites(x, y);
            EXPECT_EQ(6, nn + ns) << nn << " " << ns;
        }
    }

    const uint cx = L/2;
    const uint cy = W/2;

    int dx, dy, dz;

    //COMPLETELY FLAT

    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);

    EXPECT_EQ(0, dx); EXPECT_EQ(0, dy); EXPECT_EQ(1, dz);

    //NO NEIGHBORS
    m_solver->registerHeightChange(cx, cy, 1);
    m_solver->registerHeightChange(cx, cy, 1);

    EXPECT_EQ(5, m_solver->numberOfSurroundingSolutionSites(cx, cy));
    //0 = above => dr = (0, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //1 = left => dr = (-1, 0, 0)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ(-1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(0, dz);

    //2 = right => dr = (1, 0, 0)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(0, dz);

    //3 = bottom => dr = (0, -1, 0)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 3);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(0, dz);

    //4 = top => dr = (0, 1, 0)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 4);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 1, dy); EXPECT_EQ(0, dz);


    //ONE NEIGHBOR AT LEFT
    //increase left site to 1
    m_solver->registerHeightChange(cx-1, cy, 1);

    EXPECT_EQ(4, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    //0 = above => dr = (0, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //1 = left => dr = (-1, 0, 0) BLOCKED

    //1 = right => dr = (1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ( 1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(0, dz);

    //2 = bottom => dr = (0, -1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(0, dz);

    //3 = top => dr = (0, 1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 3);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 1, dy); EXPECT_EQ(0, dz);



    //ONE NEIGHBOR AT LEFT AND ONE AT TOP
    //increase bottom site to 2 and top site to 1
    m_solver->registerHeightChange(cx-1, cy, 1);
    m_solver->registerHeightChange(cx, cy+1, 1);

    EXPECT_EQ(3, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    //0 = above => dr = (0, 0, 2)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //1 = left => dr = (-1, 0, 1) BLOCKED

    //1 = right => dr = (1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ( 1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(0, dz);

    //2 = bottom => dr = (0, -1, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(0, dz);

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
    diffusionEvent->clearDiffusionReactions();
    //1 2 2 0 2
    EXPECT_EQ(3, m_solver->numberOfSurroundingSolutionSites(cx, cy));

    //0 = above => dr = (0, 0, 2)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 0);
    EXPECT_EQ( 0, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(1, dz);

    //1 = left => dr = (-1, 0, 1)
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 1);
    EXPECT_EQ(-1, dx); EXPECT_EQ( 0, dy); EXPECT_EQ(0, dz);

    //2 = right => dr = (1, 0, 1) //BLOCKED

    //2 = bottom => dr = (0, -1, 1) //BLOCKED
    m_solver->getSolutionSite(cx, cy, dx, dy, dz, 2);
    EXPECT_EQ( 0, dx); EXPECT_EQ(-1, dy); EXPECT_EQ(0, dz);

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
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();
    //    rng.initialize(1000);

    //    OfflatticeMonteCarloBoundary *diffusionEvent = static_cast<OfflatticeMonteCarloBoundary*>(m_diffusionEvent);

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_EQ(m_solver->nNeighbors(x, y) + m_solver->numberOfSurroundingSolutionSites(x, y, m_solver->height(x, y)), 6);
        }
    }

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
    EXPECT_EQ(2, r->calculateNumberOfFreePaths());

    //Deposit 4 paticles (neighbors of center site)
    diffusionEvent->diffusionReaction(0, 1, 2)->executeReaction(0, 0, -1);
    diffusionEvent->diffusionReaction(2, 1, 2)->executeReaction(0, 0, -1);
    diffusionEvent->diffusionReaction(1, 0, 2)->executeReaction(0, 0, -1);
    diffusionEvent->diffusionReaction(1, 2, 2)->executeReaction(0, 0, -1);
    EXPECT_EQ(6, r->calculateNumberOfFreePaths());

    EXPECT_EQ(1, m_solver->height(0, 1));
    EXPECT_EQ(1, m_solver->height(2, 1));
    EXPECT_EQ(1, m_solver->height(1, 0));
    EXPECT_EQ(1, m_solver->height(1, 2));

    //Deposit 4 particles (not neighbors of center, to neighbors of center)
    diffusionEvent->diffusionReaction(0, 0, 2)->executeReaction( 1, 0, 0);
    EXPECT_EQ(5, r->calculateNumberOfFreePaths());

    diffusionEvent->diffusionReaction(2, 2, 2)->executeReaction(-1, 0, 0);
    EXPECT_EQ(4, r->calculateNumberOfFreePaths());

    diffusionEvent->diffusionReaction(0, 2, 2)->executeReaction( 0,-1, 0);
    EXPECT_EQ(3, r->calculateNumberOfFreePaths());

    diffusionEvent->diffusionReaction(2, 0, 2)->executeReaction( 0, 1, 0);
    EXPECT_EQ(2, r->calculateNumberOfFreePaths());

    EXPECT_EQ(2, m_solver->height(0, 1));
    EXPECT_EQ(2, m_solver->height(2, 1));
    EXPECT_EQ(2, m_solver->height(1, 0));
    EXPECT_EQ(2, m_solver->height(1, 2));

    EXPECT_EQ(1, diffusionEvent->numberOfDiffusionReactions());

    r->setZ((uint)height);
    EXPECT_EQ(5, r->calculateNumberOfFreePaths());

    r->setZ(2);
    EXPECT_EQ(2, r->calculateNumberOfFreePaths());

    uint upCount = 0;
    uint downCount = 0;
    const uint N = 1000000;

    int dx, dy, dz;

    for (uint i = 0; i < N; ++i)
    {
        r->getDiffusionPath(r->calculateNumberOfFreePaths()*rng.uniform(), dx, dy, dz);

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

    EXPECT_NEAR(0.5, upProb, 1E-2);
    EXPECT_NEAR(0.5, downProb, 1E-2);

    r->executeReaction(0, 0, -1);

    EXPECT_EQ(0, diffusionEvent->numberOfDiffusionReactions());

}

TEST_F(SOSkMCTest, dissolution)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();
    //    rng.initialize(1000);

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(x, y)) << x << " " << y;
        }
    }

    if (HasFailure())
    {
        return;
    }


    EXPECT_EQ(0, diffusionEvent->numberOfDiffusionReactions());

    //dissolve a particle
    m_solver->registerHeightChange(1, 1, -1);

    //only dissolution path should be straight up
    EXPECT_EQ(1, diffusionEvent->numberOfDiffusionReactions());
    EXPECT_TRUE(nullptr!= diffusionEvent->diffusionReaction(1, 1, 1));
    diffusionEvent->clearDiffusionReactions();

    //Sites neighboring to center should have 1 paths now
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(0, 1));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(2, 1));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(1, 0));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(1, 2));

    //those who do not should still have 1
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(0, 0));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(2, 2));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(0, 2));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(2, 0));

    //dissolve again. Now we should have 2 paths for neighbors
    m_solver->registerHeightChange(1, 1, -1);
    diffusionEvent->clearDiffusionReactions();

    //Sites neighboring to center should have 2 paths now
    EXPECT_EQ(2, m_solver->numberOfSurroundingSolutionSites(0, 1));
    EXPECT_EQ(2, m_solver->numberOfSurroundingSolutionSites(2, 1));
    EXPECT_EQ(2, m_solver->numberOfSurroundingSolutionSites(1, 0));
    EXPECT_EQ(2, m_solver->numberOfSurroundingSolutionSites(1, 2));

    //those who do not should still have 1
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(0, 0));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(2, 2));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(0, 2));
    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(2, 0));


    const uint N = 1000000;

    double nSideways = 0;
    double nUp = 0;

    int dx, dy, dz;
    for (uint i = 0; i < N; ++i)
    {
        uint n = rng.uniform()*m_solver->numberOfSurroundingSolutionSites(0, 1);
        m_solver->getSolutionSite(0, 1, dx, dy, dz, n);

        EXPECT_EQ(1, abs(dx + dy + dz)) << dx << " y " << dy << " z " << dz;
        EXPECT_EQ(0, dy); EXPECT_NE(-1, dx); EXPECT_NE(-1, dz);

        if (dx == 1)
        {
            EXPECT_EQ(0, dz);

            nSideways++;
        }

        else if (dz == 1)
        {
            EXPECT_EQ(0, dx);

            nUp++;
        }

        else
        {
            EXPECT_TRUE(false);
        }

        if (HasFailure())
        {
            return;
        }

    }

    double pSideways = nSideways/N;
    double pUp = nUp/N;

    EXPECT_NEAR(0.5, pSideways, 1e-2);
    EXPECT_NEAR(0.5, pUp, 1e-2);

}


TEST_F(SOSkMCTest, SOS_discrete_interface)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();
    //    rng.initialize(1000);

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    m_solver->registerHeightChange(1,1,1);           //one peak on a flat surface
    diffusionEvent->addDiffusionReactant(1,1,3);     //a particle two lattice units away
    diffusionEvent->addDiffusionReactant(1, 1, 4);   //a particle three lattice units over peak
    diffusionEvent->addDiffusionReactant(1, 1, 5);   //a particle four lattice units over peak
    diffusionEvent->addDiffusionReactant(2, 1, 2);   //a particle next to the gap between peak and first particle

    EXPECT_EQ(4, diffusionEvent->numberOfDiffusionReactions());
    EXPECT_EQ(1, m_solver->height(1, 1));

    diffusionEvent->diffusionReaction(2, 1, 2)->executeReaction(-1, 0, 0); //move the particle into the gap.

    //this should turn all paticles into the SOS surface with a new peak of height 3
    EXPECT_EQ(0, diffusionEvent->numberOfDiffusionReactions());
    EXPECT_EQ(5, m_solver->height(1, 1));


}



TEST_F(SOSkMCTest, SOS_allSolutionSites)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 20 + rng.uniform();

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();
    //    rng.initialize(1000);

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }


    vector<vector<int>> surroundings = {{1, 1, 4},
                                        {1, 1, 6},
                                        {0, 1, 5},
                                        {2, 1, 5},
                                        {1, 0, 5},
                                        {1, 2, 5}};

    uint i = 0;

    const uint nNeighborsMax = 6;

    uvec positions(nNeighborsMax);

    //main loop over number of active particles.
    for (uint N = 1; N < nNeighborsMax; ++N)
    {

        for (uint n = 0; n < N; ++n)
        {
            positions(n) = n;
        }

        bool finished = false;
        uint lastIndex = N-1;
        while (!finished)
        {
            SOSDiffusionReaction *r = diffusionEvent->addDiffusionReactant(1, 1, 5);

            //this is the initial position after each move
            for (uint n = 0; n < N; ++n)
            {
                diffusionEvent->addDiffusionReactant(surroundings.at(positions(n)).at(0),
                                                     surroundings.at(positions(n)).at(1),
                                                     surroundings.at(positions(n)).at(2));
            }

            diffusionEvent->dump(i);
            i++;
            for (uint k = 0; k < N; ++k)
            {
                cout << positions(k) << " ";
            }
            cout << endl << "---" << endl;

            EXPECT_EQ(nNeighborsMax - N, r->calculateNumberOfFreePaths());

            int dx, dy, dz;
            for (uint p = 0; p < r->calculateNumberOfFreePaths(); ++p)
            {
                r->getDiffusionPath(p, dx, dy, dz);
                uint x = m_solver->boundaryTransform(r->x(), r->y(), r->z(), dx, 0);
                uint y = m_solver->boundaryTransform(r->x(), r->y(), r->z(), dy, 1);
                int z = r->z() + dz;

                //we expect that there are no particles along the free diffusion paths.
                EXPECT_TRUE(diffusionEvent->diffusionReaction(x, y, z) == nullptr) << "i= " <<i <<" path: " << p
                                                                                   << " dx dy dz "
                                                                                   << dx << " "
                                                                                   << dy << " "
                                                                                   << dz;
            }

            if (positions(lastIndex) == nNeighborsMax - 1)
            {

                for (int n = N-1; n >= 0; --n)
                {
                    if (n == 0)
                    {
                        finished = true;
                        break;
                    }

                    if (positions(n-1) != positions(n) - 1)
                    {
                        positions(n-1)++;
                        break;
                    }
                }

                if (!finished)
                {
                    positions(lastIndex) = positions(lastIndex-1) + 1;
                }
            }
            else
            {
                positions(lastIndex)++;
            }

            diffusionEvent->clearDiffusionReactions();

            if (HasFailure())
            {
                return;
            }
        }

    }



}




TEST_F(SOSkMCTest, SOS_diff_closeconf)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 4;
    const int iheight = (const int) height;

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {0, 0, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();
    //    rng.initialize(1000);

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    m_solver->setHeight(1, 1, iheight - 1);

    EXPECT_EQ(m_solver->numberOfSurroundingSolutionSites(1, 1), 4);

    m_solver->setHeight(1, 0, iheight - 1);
    m_solver->setHeight(1, 2, iheight - 1);
    m_solver->setHeight(0, 1, iheight - 1);
    m_solver->setHeight(2, 1, iheight - 1);

    EXPECT_EQ(m_solver->numberOfSurroundingSolutionSites(1, 1), 0);

    m_solver->setHeight(1, 0, 0);
    diffusionEvent->clearDiffusionReactions();
    m_solver->setHeight(1, 2, 0);
    diffusionEvent->clearDiffusionReactions();
    m_solver->setHeight(0, 1, 0);
    diffusionEvent->clearDiffusionReactions();
    m_solver->setHeight(2, 1, 0);
    diffusionEvent->clearDiffusionReactions();

    diffusionEvent->addDiffusionReactant(1, 0, iheight - 1);
    diffusionEvent->addDiffusionReactant(1, 2, iheight - 1);
    diffusionEvent->addDiffusionReactant(0, 1, iheight - 1);
    diffusionEvent->addDiffusionReactant(2, 1, iheight - 1);

    EXPECT_EQ(m_solver->numberOfSurroundingSolutionSites(1, 1), 0);

}




TEST_F(SOSkMCTest, SOS_diff_mixedboundary)
{
    const uint L = 3;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 4;
    const int iheight = (const int) height;

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {2, 1, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();
    rng.initialize(1001230);

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    m_solver->setHeight(0, 2, iheight - 1);
    diffusionEvent->addDiffusionReactant(0, 1, iheight - 1);
    diffusionEvent->addDiffusionReactant(0, 0, iheight - 1);

    diffusionEvent->Diffusion::dump(0);
    diffusionEvent->dump(0);

    EXPECT_EQ(2, m_solver->numberOfSurroundingSolutionSites(0, 2));



}



TEST_F(SOSkMCTest, SOS_diff_closewall_trapped)
{
    const uint L = 10;
    const uint W = 3;
    const double alpha = 1.0;
    const double mu = 0;
    const double height = 7;

    m_solver = new SOSSolver(L, W, alpha, mu);
    setBoundariesFromIDs(m_solver, {2, 1, 0, 0}, L, W);
    m_pressureWallEvent = new FixedSurface(*m_solver, height);
    LatticeDiffusion *diffusionEvent = new LatticeDiffusion(*m_solver);
    m_diffusionEvent = diffusionEvent;

    SetUp_yo();

    primeSolver();

    diffusionEvent->clearDiffusionReactions();

    for (uint x = 0; x < L; ++x)
    {
        for (uint y = 0; y < W; ++y)
        {
            m_solver->setHeight(x, y, 0);
            diffusionEvent->clearDiffusionReactions();
        }
    }

    m_solver->setHeight(5, 1, 6);
    m_solver->setHeight(6, 1, 5);
    m_solver->setHeight(5, 2, 5);

    diffusionEvent->addDiffusionReactant(5, 0, 6);
    diffusionEvent->addDiffusionReactant(4, 1, 5);
    diffusionEvent->addDiffusionReactant(5, 2, 6);

    EXPECT_EQ(1, m_solver->numberOfSurroundingSolutionSites(5, 1));
    EXPECT_FALSE(m_solver->surfaceReaction(5, 1).rateExpression() == 0);

}









