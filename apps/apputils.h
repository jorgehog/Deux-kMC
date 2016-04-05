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

class SurfaceCounter : public Event<uint>, public Observer<Subjects>
{
public:
    SurfaceCounter(const uint totalCount) :
        Event("SC", "%", true),
        Observer(),
        m_totalCount(totalCount)
    {

    }

private:

    const uint m_totalCount;
    uint m_counter;

    // Observer interface
public:
    void initializeObserver(const Subjects &subject)
    {
        (void) subject;
        m_counter = 0;
    }

    void notifyObserver(const Subjects &subject)
    {
        (void) subject;
        m_counter++;
    }

    // Event interface
public:
    void execute()
    {
        setValue(m_counter/double(m_totalCount)*100);

        if (m_counter >= m_totalCount)
        {
            terminateLoop("Counter exceeds limit");
        }
    }
};


void initializeSurface(SOSSolver &solver, const string type,
                       const uint diffusionInt = 0,
                       const uint surfaceThermCycles = 10000,
                       const bool output = true)
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

        setBoundariesFromIDs(&thermSolver, {0,0,0,0}, L, W);

        SurfaceCounter counter(surfaceThermCycles);
        thermSolver.registerObserver(&counter);

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

        lattice.addEvent(thermSolver);
        lattice.addEvent(conf);
        lattice.addEvent(counter);

        if (conc->hasDiscreteParticles())
        {
            lattice.addEvent(pnc);
            conc->setDependency(pnc);
        }

        lattice.addEvent(conc);

        initializeSurface(thermSolver, "random");

        const uint every = surfaceThermCycles/10;
        lattice.enableOutput(output, every);
        lattice.enableProgressReport(false);
        lattice.enableEventValueStorage(false, false);


        lattice.eventLoop(1000*surfaceThermCycles);

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


inline double getMFPTConstant(const double h0, const double alpha, const int type)
{
    using maptype = std::map<pair<uint, uint>, double>;

    const maptype valueMapRadial = {{{ 0,  0}, 0.080},
                                    {{ 0,  1}, 0.078},
                                    {{ 0,  2}, 0.074},
                                    {{ 0,  3}, 0.070},
                                    {{ 0,  4}, 0.066},
                                    {{ 0,  5}, 0.065},
                                    {{ 0,  6}, 0.068},
                                    {{ 0,  7}, 0.066},
                                    {{ 0,  8}, 0.066},
                                    {{ 0,  9}, 0.063},
                                    {{ 0, 10}, 0.071},
                                    {{ 0, 11}, 0.069},
                                    {{ 0, 12}, 0.074},
                                    {{ 0, 13}, 0.073},
                                    {{ 0, 14}, 0.081},
                                    {{ 0, 15}, 0.119},
                                    {{ 1,  0}, 0.077},
                                    {{ 1,  1}, 0.073},
                                    {{ 1,  2}, 0.072},
                                    {{ 1,  3}, 0.069},
                                    {{ 1,  4}, 0.069},
                                    {{ 1,  5}, 0.066},
                                    {{ 1,  6}, 0.066},
                                    {{ 1,  7}, 0.068},
                                    {{ 1,  8}, 0.064},
                                    {{ 1,  9}, 0.072},
                                    {{ 1, 10}, 0.069},
                                    {{ 1, 11}, 0.076},
                                    {{ 1, 12}, 0.078},
                                    {{ 1, 13}, 0.073},
                                    {{ 1, 14}, 0.085},
                                    {{ 1, 15}, 0.110},
                                    {{ 2,  0}, 0.078},
                                    {{ 2,  1}, 0.078},
                                    {{ 2,  2}, 0.072},
                                    {{ 2,  3}, 0.069},
                                    {{ 2,  4}, 0.068},
                                    {{ 2,  5}, 0.069},
                                    {{ 2,  6}, 0.067},
                                    {{ 2,  7}, 0.073},
                                    {{ 2,  8}, 0.071},
                                    {{ 2,  9}, 0.068},
                                    {{ 2, 10}, 0.067},
                                    {{ 2, 11}, 0.072},
                                    {{ 2, 12}, 0.080},
                                    {{ 2, 13}, 0.084},
                                    {{ 2, 14}, 0.092},
                                    {{ 2, 15}, 0.105},
                                    {{ 3,  0}, 0.093},
                                    {{ 3,  1}, 0.077},
                                    {{ 3,  2}, 0.073},
                                    {{ 3,  3}, 0.081},
                                    {{ 3,  4}, 0.075},
                                    {{ 3,  5}, 0.076},
                                    {{ 3,  6}, 0.071},
                                    {{ 3,  7}, 0.075},
                                    {{ 3,  8}, 0.073},
                                    {{ 3,  9}, 0.073},
                                    {{ 3, 10}, 0.080},
                                    {{ 3, 11}, 0.077},
                                    {{ 3, 12}, 0.096},
                                    {{ 3, 13}, 0.097},
                                    {{ 3, 14}, 0.103},
                                    {{ 3, 15}, 0.105}};

    const maptype valueMapPathfinding = {{{ 0,  0}, 0.089},
                                         {{ 0,  1}, 0.078},
                                         {{ 0,  2}, 0.078},
                                         {{ 0,  3}, 0.072},
                                         {{ 0,  4}, 0.068},
                                         {{ 0,  5}, 0.073},
                                         {{ 0,  6}, 0.068},
                                         {{ 0,  7}, 0.071},
                                         {{ 0,  8}, 0.069},
                                         {{ 0,  9}, 0.079},
                                         {{ 0, 10}, 0.076},
                                         {{ 0, 11}, 0.071},
                                         {{ 0, 12}, 0.073},
                                         {{ 0, 13}, 0.084},
                                         {{ 0, 14}, 0.097},
                                         {{ 0, 15}, 0.101},
                                         {{ 1,  0}, 0.085},
                                         {{ 1,  1}, 0.080},
                                         {{ 1,  2}, 0.077},
                                         {{ 1,  3}, 0.073},
                                         {{ 1,  4}, 0.080},
                                         {{ 1,  5}, 0.071},
                                         {{ 1,  6}, 0.070},
                                         {{ 1,  7}, 0.071},
                                         {{ 1,  8}, 0.076},
                                         {{ 1,  9}, 0.078},
                                         {{ 1, 10}, 0.079},
                                         {{ 1, 11}, 0.085},
                                         {{ 1, 12}, 0.080},
                                         {{ 1, 13}, 0.087},
                                         {{ 1, 14}, 0.099},
                                         {{ 1, 15}, 0.114},
                                         {{ 2,  0}, 0.078},
                                         {{ 2,  1}, 0.086},
                                         {{ 2,  2}, 0.080},
                                         {{ 2,  3}, 0.081},
                                         {{ 2,  4}, 0.078},
                                         {{ 2,  5}, 0.084},
                                         {{ 2,  6}, 0.071},
                                         {{ 2,  7}, 0.071},
                                         {{ 2,  8}, 0.085},
                                         {{ 2,  9}, 0.085},
                                         {{ 2, 10}, 0.080},
                                         {{ 2, 11}, 0.085},
                                         {{ 2, 12}, 0.091},
                                         {{ 2, 13}, 0.092},
                                         {{ 2, 14}, 0.099},
                                         {{ 2, 15}, 0.123},
                                         {{ 3,  0}, 0.089},
                                         {{ 3,  1}, 0.087},
                                         {{ 3,  2}, 0.088},
                                         {{ 3,  3}, 0.087},
                                         {{ 3,  4}, 0.078},
                                         {{ 3,  5}, 0.081},
                                         {{ 3,  6}, 0.078},
                                         {{ 3,  7}, 0.080},
                                         {{ 3,  8}, 0.082},
                                         {{ 3,  9}, 0.080},
                                         {{ 3, 10}, 0.089},
                                         {{ 3, 11}, 0.092},
                                         {{ 3, 12}, 0.103},
                                         {{ 3, 13}, 0.114},
                                         {{ 3, 14}, 0.116},
                                         {{ 3, 15}, 0.141}};

    const maptype* valueMap;
    if (type == 0)
    {
        valueMap = &valueMapRadial;
    }

    else
    {
        valueMap = &valueMapPathfinding;
    }

    const vec knownHeights = {5, 10, 15, 20};

    const vec knownAlphas = {0.5, 0.6, 0.7, 0.8,
                             0.9, 1, 1.1, 1.2,
                             1.3, 1.4, 1.5, 1.6,
                             1.7, 1.8, 1.9, 2};

    uint ia = 0;
    uint ih = 0;

    while (h0 != knownHeights(ih))
    {
        ih++;

        if (ih == knownHeights.size())
        {
            throw std::runtime_error("h0 not found.");
        }
    }

    while (alpha > knownAlphas(ia))
    {
        ia++;
    }

    double c0, c1;
    double da;
    double daa;

    if (valueMap->find({ih, ia}) == valueMap->end())
    {
        throw std::runtime_error("combination not found.");
    }

    if (ia != 0)
    {
        c0 = valueMap->at({ih, ia-1});
        c1 = valueMap->at({ih, ia});

        da = alpha - knownAlphas(ia-1);
        daa = knownAlphas(ia) - knownAlphas(ia-1);
    }
    else
    {
        c0 = valueMap->at({ih, ia});
        c1 = valueMap->at({ih, ia+1});

        da = alpha - knownAlphas(ia);
        daa = knownAlphas(ia+1) - knownAlphas(ia);
    }

    return c0 + (c1-c0)*da/daa;

}
























