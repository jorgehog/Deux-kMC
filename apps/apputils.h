#pragma once

#include <SOSkMC.h>

#include <sstream>
#include <BADAss/badass.h>
#include <ignis/include/ignis.h>

using namespace std;

uint getProc(int argv, char** argc)
{
    if (argv == 1)
    {
        return 0u;
    }

    else
    {
        return uint(atoi(argc[1]));
    }
}

string getTail(int argv, char** argc)
{
    string ending;

    if (argv == 1)
    {
        ending = "";
    }

    else
    {
        int proc = getProc(argv, argc);

        stringstream procEnding;
        procEnding << "_" << proc;

        ending = procEnding.str();
    }

    return ending;

}

string addProcEnding(int argv, char** argc, string filename, string ending)
{
    stringstream s;
    s << filename << getTail(argv, argc) << "." << ending;

    return s.str();
}

string getCfgName(int argv, char** argc)
{
    string cfgName;

    if (argv != 1)
    {
        BADAss(argv, ==, 3, "Usage: executable proc cfgname");

        cfgName = argc[2];
    }

    else
    {
        cfgName = "infiles/spiralgrowth.cfg";
    }

    return cfgName;
}

class EquilibrationOrganizer
{
public:
    EquilibrationOrganizer(MainMesh<uint> &lattice,
                           const bool equilibrate,
                           const bool reset,
                           const bool skipEvents) :
        m_lattice(lattice),
        m_equilibrate(equilibrate),
        m_reset(reset),
        m_skipEvents(skipEvents)
    {

    }

    ~EquilibrationOrganizer()
    {
        m_neglectedEvents.clear();
    }

    void prepare(vector<Event<uint>*> equilibriationEvents,
                 vector<Event<uint>*> staticEvents)
    {
        m_equilibrationEvents = equilibriationEvents;

        if (m_equilibrate && m_skipEvents)
        {

            //Making sure we do not remove any events that the equlibration
            //events might depend on
            vector<uint> dependencies;
            for (uint i = 0; i < m_lattice.getEvents().size(); ++i)
            {
                const Event<uint> *event = m_lattice.getEvents().at(i);

                //If the event is not to be removed.
                if (std::find(staticEvents.begin(), staticEvents.end(), event) != staticEvents.end())
                {
                    dependencies.push_back(i);
                    continue;
                }

                for (const Event<uint> *staticEvent : staticEvents)
                {
                    if (staticEvent->dependsOn(event))
                    {
                        dependencies.push_back(i);
                        break;
                    }
                }

                for (const Event<uint> *eqEvent : m_equilibrationEvents)
                {
                    if (eqEvent->dependsOn(event))
                    {
                        dependencies.push_back(i);
                        break;
                    }
                }
            }

            //Keep the events which we do not use
            for (uint i = 0; i < m_lattice.getEvents().size(); ++i)
            {
                if (std::find(dependencies.begin(), dependencies.end(), i) != dependencies.end())
                {
                    continue;
                }

                m_neglectedEvents.push_back(m_lattice.getEvents().at(i));
            }

            //remove events for event lists
            for (int i = m_lattice.getEvents().size() - 1; i >= 0; --i)
            {
                if (std::find(dependencies.begin(), dependencies.end(), i) != dependencies.end())
                {
                    continue;
                }

                m_lattice.removeEvent(i);
            }

            //Add the equilibration events
            for (Event<uint> *eqEvent : m_equilibrationEvents)
            {
                m_lattice.addEvent(eqEvent);
            }
        }
    }

    void reset(const uint nCycles)
    {
        if (m_equilibrate && m_reset)
        {
            //Remove the equilibriation events
            for (Event<uint> *eqEvent : m_equilibrationEvents)
            {
                m_lattice.removeEvent(eqEvent);
            }

            //add the skipped events back in
            if (m_skipEvents)
            {
                for (auto &event : m_neglectedEvents)
                {
                    m_lattice.addEvent(event);
                }
            }

            //restart the simulation
            m_lattice.eventLoop(nCycles);
        }
    }

private:

    MainMesh<uint> &m_lattice;
    const bool m_equilibrate;
    const bool m_reset;
    const bool m_skipEvents;

    vector<LatticeEvent*> m_neglectedEvents;
    vector<LatticeEvent*> m_equilibrationEvents;

};

const Boundary *getBoundaryFromID(SOSSolver *solver,
                                  const uint ID,
                                  const uint dim,
                                  const uint span,
                                  const uint yspan,
                                  Boundary::orientations orientation,
                                  const int boundaryHeight = 0,
                                  const uint averageHeightDepth = 0)
{
    uint location = orientation == Boundary::orientations::FIRST ? 0 : (span - 1);

    switch (ID) {
    case 0:
        return new Periodic(span, orientation);
        break;
    case 1:
        return new Edge(location, orientation);
        break;
    case 2:
        return new Open(orientation);
        break;
    case 3:
        return new Reflecting(location, orientation);
        break;
    case 4:
        return new ConstantHeight(boundaryHeight, location, orientation);
    case 5:
        return new AverageHeightBoundary(*solver, averageHeightDepth, dim, span, yspan, orientation, location);
    default:
        cerr << "invalid boundary: " << ID << endl;
        return nullptr;
        break;
    }

}

void setBoundariesFromIDs(SOSSolver *solver,
                          const vector<uint> IDs,
                          const uint L, const uint W,
                          const int boundaryHeight = 0,
                          const uint averageHeightDepth = 0)
{
    const Boundary* leftBoundary = getBoundaryFromID(solver, IDs.at(0), 0, L, W, Boundary::orientations::FIRST, boundaryHeight, averageHeightDepth);
    const Boundary* rightBoundary = getBoundaryFromID(solver, IDs.at(1), 0, L, W, Boundary::orientations::LAST, boundaryHeight, averageHeightDepth);
    const Boundary* bottomBoundary = getBoundaryFromID(solver, IDs.at(2), 1, W, L, Boundary::orientations::FIRST, boundaryHeight, averageHeightDepth);
    const Boundary* topBoundary = getBoundaryFromID(solver, IDs.at(3), 1, W, L, Boundary::orientations::LAST, boundaryHeight, averageHeightDepth);

    solver->setBoundaries({{leftBoundary, rightBoundary},
                           {bottomBoundary, topBoundary}});
}



























