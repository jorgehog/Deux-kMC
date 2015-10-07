#pragma once

#include <SOSkMC.h>
#include <gtest/gtest.h>
#include <sys/time.h>

#include "../apputils.h"

#define deletehaxx(name) if (name != nullptr) delete name; name = nullptr

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
        rng.initialize(time(nullptr));
    }

    void SetUp_yo()
    {

        // Code here will be called immediately after the constructor (right
        // before each test).


        m_lattice = new Lattice();

        m_lattice->addEvent(m_solver);
        m_lattice->addEvent(m_pressureWallEvent);
        m_lattice->addEvent(m_diffusionEvent);

        m_averageHeight = new AverageHeight(*m_solver);
        m_lattice->addEvent(m_averageHeight);

        m_lattice->enableOutput(false);
        m_lattice->enableEventValueStorage(false, false);

    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test (right
        // before the destructor).

        deletehaxx(m_pressureWallEvent);
        deletehaxx(m_diffusionEvent);
        deletehaxx(m_averageHeight);
        deletehaxx(m_solver);
        deletehaxx(m_lattice);

    }

    // Objects declared here can be used by all tests in the test case for Foo.

    void primeSolver()
    {
        BasicInitializeEvent<uint> primer("kmc primer", [&] (BasicInitializeEvent<uint> *event)
        {
                event->stopLoop();
        });

        m_lattice->addEvent(primer);

        m_lattice->eventLoop(1);

        m_lattice->removeEvent(&primer);
    }

    SOSSolver &solver()
    {
        return *m_solver;
    }

    SOSSolver *m_solver = nullptr;
    AverageHeight *m_averageHeight = nullptr;
    ConfiningSurface *m_pressureWallEvent = nullptr;
    Diffusion *m_diffusionEvent = nullptr;

    Lattice *m_lattice = nullptr;

};
