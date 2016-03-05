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

string getCfgName(int argv, char** argc, string name = "spiralgrowth")
{
    string cfgName;

    if (argv != 1)
    {
        BADAss(argv, ==, 3, "Usage: executable proc cfgname");

        cfgName = argc[2];
    }

    else
    {
        stringstream ss;
        ss << "infiles/" << name << ".cfg";
        cfgName = ss.str();
    }

    return cfgName;
}


template<typename eT>
class EventIsolator
{
    typedef Event<eT>* EeT;

public:
    EventIsolator(MeshField<eT> &lattice) :
        m_lattice(lattice),
        m_skippedEvents([] (const EeT a, const EeT b) {return a->meshAddress() < b->meshAddress();})
    {

    }

    void isolate(vector<string> eventNames)
    {
        vector<Event<eT>*> events;

        for (Event<eT> *event : m_lattice.getEvents())
        {
            if (std::find(eventNames.begin(), eventNames.end(), event->type()) != eventNames.end())
            {
                events.push_back(event);
            }
        }

        isolate(events);
    }

    void isolate(vector<Event<eT>*> events)
    {

        //events are added in the order in which they appear in the mainmesh.

        set<uint> isolatedEvents;

        for (uint i = 0; i < events.size(); ++i)
        {
            Event<eT> *isolatedEvent = events.at(i);

            for (uint j = 0; j < m_lattice.getEvents().size(); ++j)
            {
                Event<eT> *event = m_lattice.getEvents().at(j);

                //if the event is not in the isolated event list.
                if (std::find(events.begin(), events.end(), event) == events.end())
                {

                    //and the isolated event depends on it
                    if (isolatedEvent->dependsOn(event))
                    {
                        isolatedEvents.insert(j);
                    }

                    else
                    {
                        m_skippedEvents.insert(event);
                    }
                }

                //else we append it
                else
                {
                    isolatedEvents.insert(j);
                }
            }
        }

        //remove events for event lists
        for (int i = m_lattice.getEvents().size() - 1; i >= 0; --i)
        {
            uint I = i;

            if (std::find(isolatedEvents.begin(), isolatedEvents.end(), I) != isolatedEvents.end())
            {
                continue;
            }

            m_lattice.removeEvent(i);
        }

    }

    void release()
    {
        for (Event<eT> *skippedEvent : m_skippedEvents)
        {
            m_lattice.addEvent(skippedEvent);
        }

        m_skippedEvents.clear();
    }

private:

    MeshField<eT> &m_lattice;
    set<EeT, std::function<bool(EeT, EeT)>> m_skippedEvents;
};


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

void initializeSurface(SOSSolver &solver, const string type,
                       const uint diffusionInt = 0,
                       const uint surfaceThermCycles = 10000)
{
    int maxSpan;

    if (solver.confiningSurfaceEvent().hasSurface())
    {
        maxSpan = solver.confiningSurfaceEvent().height()/10;

        if (maxSpan == 0)
        {
            maxSpan = 1;
        }
    }
    else
    {
        maxSpan = 3;
    }

    const uint L = solver.length();
    const uint W = solver.width();

    if (type == "fracture")
    {
        int noise;
        double value = 0;

        int TwoThirdsWay = (2*solver.confiningSurfaceEvent().height())/3;

        double delta = TwoThirdsWay/double(solver.length());

        for (uint x = 0; x < solver.length(); ++x)
        {
            value = value + delta;
            for (uint y = 0; y < solver.width(); ++y)
            {
                noise = -1 + 3*rng.uniform();
                solver.setHeight(x, y, value + noise, false);
            }
        }
    }

    else if (type == "flat")
    {
        solver.setHeights(zeros<imat>(solver.length(), solver.width()), false);
    }

    else if (type == "singlestep")
    {
        for (uint x = 0; x < L; ++x)
        {
            for (uint y = 0; y < W; ++y)
            {
                if (x < L/10 ||
                        x > (9*L)/10 ||
                        y < W/10 ||
                        y > (9*W)/10)
                {
                    solver.setHeight(x, y, 1, false);
                }

                else
                {
                    solver.setHeight(x, y, 0, false);
                }
            }
        }
    }

    else if (type == "random")
    {

        solver.setHeights(randi(solver.length(),
                                solver.width(),
                                distr_param(-maxSpan, maxSpan)),
                          false);

        int m = accu(solver.heights());

        uint x;
        uint y;
        if (m != 0)
        {
            int N = abs(m);
            int direction = -m/N;

            for (int shift = 0; shift < N; ++shift)
            {
                do
                {
                    x = rng.uniform()*solver.length();
                    y = rng.uniform()*solver.width();
                } while (abs(solver.height(x, y)) == maxSpan);

                solver.setHeight(x, y, solver.height(x, y) + direction, false);
            }
        }

    }

    else if (type == "temple")
    {
        const uint templeSize = 2*maxSpan + 1;

        for (uint x = 0; x < L; ++x)
        {
            const int templeIndexX = x%templeSize;
            const int xrank = maxSpan - abs(maxSpan - templeIndexX);

            const int polarityX = (x/templeSize)%2 == 0 ? 1 : -1;

            for (uint y = 0; y < W; ++y)
            {
                const int templeIndexY = y%templeSize;
                const int yrank = maxSpan - abs(maxSpan - templeIndexY);

                const int polarityY = (y/templeSize)%2 == 0 ? 1 : -1;

                const int h = xrank < yrank ? xrank : yrank;

                const int polarity = polarityX*polarityY;

                const int shift = polarity < 0 ? -1 : 0;

                solver.setHeight(x, y, shift + polarity*h);
            }
        }

    }

    else if (type == "thermalized")
    {
        Lattice lattice;
        SOSSolver thermSolver(L, W, solver.alpha(), solver.gamma(), solver.surfaceDiffusion());
        Diffusion *conc;
        ConfiningSurface *conf;

        if (solver.confiningSurfaceEvent().hasSurface())
        {
            conf = new ConstantConfinement(thermSolver, solver.confiningSurfaceEvent().height());
        }

        else
        {
            conf = new NoConfinement(thermSolver);
        }

        bool onlattice = diffusionInt == 2;

        if (onlattice)
        {
            conc = new LatticeDiffusion(thermSolver);
        }

        else if (diffusionInt == 3)
        {
            RadialFirstPassage *fpce = dynamic_cast<RadialFirstPassage*>(&solver.diffusionEvent());
            conc = new RadialFirstPassage(thermSolver, fpce->maxdt(), fpce->depositionBoxHalfSize(), fpce->c());
        }

        else if (diffusionInt == 4)
        {
            AStarFirstPassage *fpce = dynamic_cast<AStarFirstPassage*>(&solver.diffusionEvent());
            conc = new AStarFirstPassage(thermSolver, fpce->maxdt(), fpce->depositionBoxHalfSize(), fpce->c());
        }

        else
        {
            conc = new ConstantConcentration(thermSolver);
        }

        ParticleNumberConservator pnc(thermSolver);

        conf->registerObserver(conc);

        setBoundariesFromIDs(&thermSolver, {0,0,0,0}, L, W);

        lattice.addEvent(thermSolver);
        lattice.addEvent(conf);

        if (conc->hasDiscreteParticles())
        {
            lattice.addEvent(pnc);
        }

        lattice.addEvent(conc);

        initializeSurface(thermSolver, "random");

        const uint every = surfaceThermCycles/10;
        lattice.enableOutput(true, every);
        lattice.enableProgressReport();
        lattice.enableEventValueStorage(false, false);

        lattice.eventLoop(surfaceThermCycles);

        solver.setHeights(thermSolver.heights(), false);

        if (solver.confiningSurfaceEvent().hasSurface())
        {
            const double &hPrev = solver.confiningSurfaceEvent().height();
            const double hMean = thermSolver.averageHeight();

            solver.confiningSurfaceEvent().setHeight(hPrev + hMean);
        }

        if (onlattice)
        {
            LatticeDiffusion* diffEvent = dynamic_cast<LatticeDiffusion*>(&solver.diffusionEvent());
            for (const auto &m : dynamic_cast<LatticeDiffusion*>(conc)->diffusionReactionsMap())
            {
                SOSDiffusionReaction *r = m.second;
                diffEvent->addDiffusionReactant(r->x(), r->y(), r->z());
            }
        }

        else if (diffusionInt == 4 || diffusionInt == 3)
        {
            OfflatticeMonteCarlo *solverOfflattice = dynamic_cast<OfflatticeMonteCarlo*>(&solver.diffusionEvent());
            OfflatticeMonteCarlo *thermOfflattice = dynamic_cast<OfflatticeMonteCarlo*>(conc);

            for (uint n = 0; n < thermOfflattice->numberOfParticles(); ++n)
            {
                solverOfflattice->insertParticle(thermOfflattice->particlePositions(0, n),
                                                 thermOfflattice->particlePositions(1, n),
                                                 thermOfflattice->particlePositions(2, n));
            }

            BADAssEqual(solverOfflattice->numberOfParticles(), thermOfflattice->numberOfParticles());
        }

        delete conc;
        delete conf;
    }

    else if (type == "none")
    {
        return;
    }

    else
    {
        throw std::runtime_error("invalid surfacetype: " + type);
    }



}

























