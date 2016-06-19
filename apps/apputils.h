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

Boundary *getBoundaryFromID(SOSSolver *solver,
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
    case 6:
        return new ReflectingSurfaceOpenSolution(location, orientation);
    default:
        cerr << "invalid boundary: " << ID << endl;
        return nullptr;
        break;
    }

}

vector<vector<Boundary *> >
setBoundariesFromIDs(SOSSolver *solver,
                     const vector<uint> IDs,
                     const uint L, const uint W,
                     const int boundaryHeight = 0,
                     const uint averageHeightDepth = 0)
{
    Boundary* leftBoundary = getBoundaryFromID(solver, IDs.at(0), 0, L, W, Boundary::orientations::FIRST, boundaryHeight, averageHeightDepth);
    Boundary* rightBoundary = getBoundaryFromID(solver, IDs.at(1), 0, L, W, Boundary::orientations::LAST, boundaryHeight, averageHeightDepth);
    Boundary* bottomBoundary = getBoundaryFromID(solver, IDs.at(2), 1, W, L, Boundary::orientations::FIRST, boundaryHeight, averageHeightDepth);
    Boundary* topBoundary = getBoundaryFromID(solver, IDs.at(3), 1, W, L, Boundary::orientations::LAST, boundaryHeight, averageHeightDepth);

    vector<vector<Boundary *> > boundaries = {{leftBoundary, rightBoundary},
                                              {bottomBoundary, topBoundary}};

    vector<vector<const Boundary *> > boundaries2 = {{leftBoundary, rightBoundary},
                                                     {bottomBoundary, topBoundary}};

    solver->setBoundaries(boundaries2);

    return boundaries;
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

        vector<vector<Boundary *> > b = setBoundariesFromIDs(&thermSolver, {0,0,0,0}, L, W);

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

        delete b[0][0];
        delete b[0][1];
        delete b[1][0];
        delete b[1][1];
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

    const maptype valueMapRadial = {{{ 0,  0}, 0.048},
                                    {{ 0,  1}, 0.047},
                                    {{ 0,  2}, 0.042},
                                    {{ 0,  3}, 0.041},
                                    {{ 0,  4}, 0.041},
                                    {{ 0,  5}, 0.041},
                                    {{ 0,  6}, 0.045},
                                    {{ 0,  7}, 0.041},
                                    {{ 0,  8}, 0.043},
                                    {{ 0,  9}, 0.044},
                                    {{ 0, 10}, 0.044},
                                    {{ 0, 11}, 0.046},
                                    {{ 0, 12}, 0.047},
                                    {{ 0, 13}, 0.050},
                                    {{ 0, 14}, 0.054},
                                    {{ 0, 15}, 0.059},
                                    {{ 1,  0}, 0.064},
                                    {{ 1,  1}, 0.063},
                                    {{ 1,  2}, 0.057},
                                    {{ 1,  3}, 0.057},
                                    {{ 1,  4}, 0.057},
                                    {{ 1,  5}, 0.056},
                                    {{ 1,  6}, 0.058},
                                    {{ 1,  7}, 0.059},
                                    {{ 1,  8}, 0.058},
                                    {{ 1,  9}, 0.055},
                                    {{ 1, 10}, 0.057},
                                    {{ 1, 11}, 0.063},
                                    {{ 1, 12}, 0.066},
                                    {{ 1, 13}, 0.062},
                                    {{ 1, 14}, 0.082},
                                    {{ 1, 15}, 0.082},
                                    {{ 2,  0}, 0.071},
                                    {{ 2,  1}, 0.068},
                                    {{ 2,  2}, 0.066},
                                    {{ 2,  3}, 0.064},
                                    {{ 2,  4}, 0.068},
                                    {{ 2,  5}, 0.063},
                                    {{ 2,  6}, 0.057},
                                    {{ 2,  7}, 0.063},
                                    {{ 2,  8}, 0.063},
                                    {{ 2,  9}, 0.068},
                                    {{ 2, 10}, 0.071},
                                    {{ 2, 11}, 0.064},
                                    {{ 2, 12}, 0.078},
                                    {{ 2, 13}, 0.078},
                                    {{ 2, 14}, 0.090},
                                    {{ 2, 15}, 0.085},
                                    {{ 3,  0}, 0.081},
                                    {{ 3,  1}, 0.069},
                                    {{ 3,  2}, 0.073},
                                    {{ 3,  3}, 0.071},
                                    {{ 3,  4}, 0.069},
                                    {{ 3,  5}, 0.068},
                                    {{ 3,  6}, 0.066},
                                    {{ 3,  7}, 0.072},
                                    {{ 3,  8}, 0.071},
                                    {{ 3,  9}, 0.072},
                                    {{ 3, 10}, 0.082},
                                    {{ 3, 11}, 0.076},
                                    {{ 3, 12}, 0.079},
                                    {{ 3, 13}, 0.092},
                                    {{ 3, 14}, 0.093},
                                    {{ 3, 15}, 0.109}};

    const maptype valueMapPathfind = {{{ 0,  0}, 0.057},
                                      {{ 0,  1}, 0.051},
                                      {{ 0,  2}, 0.047},
                                      {{ 0,  3}, 0.046},
                                      {{ 0,  4}, 0.044},
                                      {{ 0,  5}, 0.044},
                                      {{ 0,  6}, 0.044},
                                      {{ 0,  7}, 0.044},
                                      {{ 0,  8}, 0.046},
                                      {{ 0,  9}, 0.045},
                                      {{ 0, 10}, 0.045},
                                      {{ 0, 11}, 0.047},
                                      {{ 0, 12}, 0.052},
                                      {{ 0, 13}, 0.052},
                                      {{ 0, 14}, 0.061},
                                      {{ 0, 15}, 0.064},
                                      {{ 1,  0}, 0.079},
                                      {{ 1,  1}, 0.071},
                                      {{ 1,  2}, 0.061},
                                      {{ 1,  3}, 0.062},
                                      {{ 1,  4}, 0.057},
                                      {{ 1,  5}, 0.057},
                                      {{ 1,  6}, 0.062},
                                      {{ 1,  7}, 0.062},
                                      {{ 1,  8}, 0.066},
                                      {{ 1,  9}, 0.060},
                                      {{ 1, 10}, 0.062},
                                      {{ 1, 11}, 0.076},
                                      {{ 1, 12}, 0.068},
                                      {{ 1, 13}, 0.082},
                                      {{ 1, 14}, 0.083},
                                      {{ 1, 15}, 0.085},
                                      {{ 2,  0}, 0.072},
                                      {{ 2,  1}, 0.073},
                                      {{ 2,  2}, 0.071},
                                      {{ 2,  3}, 0.064},
                                      {{ 2,  4}, 0.064},
                                      {{ 2,  5}, 0.066},
                                      {{ 2,  6}, 0.068},
                                      {{ 2,  7}, 0.068},
                                      {{ 2,  8}, 0.071},
                                      {{ 2,  9}, 0.072},
                                      {{ 2, 10}, 0.068},
                                      {{ 2, 11}, 0.079},
                                      {{ 2, 12}, 0.082},
                                      {{ 2, 13}, 0.099},
                                      {{ 2, 14}, 0.096},
                                      {{ 2, 15}, 0.102},
                                      {{ 3,  0}, 0.076},
                                      {{ 3,  1}, 0.083},
                                      {{ 3,  2}, 0.081},
                                      {{ 3,  3}, 0.073},
                                      {{ 3,  4}, 0.070},
                                      {{ 3,  5}, 0.073},
                                      {{ 3,  6}, 0.082},
                                      {{ 3,  7}, 0.074},
                                      {{ 3,  8}, 0.082},
                                      {{ 3,  9}, 0.076},
                                      {{ 3, 10}, 0.082},
                                      {{ 3, 11}, 0.079},
                                      {{ 3, 12}, 0.092},
                                      {{ 3, 13}, 0.096},
                                      {{ 3, 14}, 0.127},
                                      {{ 3, 15}, 0.127}};

    const vec knownHeights = {5, 10, 15, 20};

    const vec knownAlphas = {0.5, 0.6, 0.7, 0.8,
                             0.9, 1, 1.1, 1.2,
                             1.3, 1.4, 1.5, 1.6,
                             1.7, 1.8, 1.9, 2};



    const maptype* valueMap;
    if (type == 0)
    {
        valueMap = &valueMapRadial;
    }

    else
    {
        valueMap = &valueMapPathfind;
    }

    const uvec knownAlpha = arma::find(knownAlphas == alpha);
    const uvec knownHeight = arma::find(knownHeights == h0);

    const bool alphaKnown = knownAlpha.size() != 0;
    const bool hKnown = knownHeight.size() != 0;

    uint ia = 0;

    if (alphaKnown)
    {
        ia = knownAlpha(0);
    }
    else
    {
        return 0.1;
    }

    uint ih = 0;

    if (hKnown)
    {
        ih = knownHeight(0);
        return valueMap->at({ih, ia});
    }
    else
    {
        return 0.1;
    }

}

H5Wrapper::Member &setuph5(H5Wrapper::Root &h5root, const uint proc, const uint L, const uint W)
{
    stringstream sizeDesc;
    sizeDesc << L << "x" << W;
    H5Wrapper::Member &sizeRoot = h5root.addMember(sizeDesc.str());

    timeval tv;
    gettimeofday(&tv, nullptr);

    __int64_t run_ID = 1000*tv.tv_sec + tv.tv_usec/1000 + 10000000000000u*proc;

    return sizeRoot.addMember(run_ID);
}

class StoreHeights : public SOSEvent
{
public:
    StoreHeights(const SOSSolver &solver, const uint interval, H5Wrapper::Member &h5group) :
        SOSEvent(solver, "storeheights"),
        m_interval(interval),
        m_h5group(h5group.addMember("stored_heights"))
    {

    }

private:
    const uint m_interval;
    H5Wrapper::Member &m_h5group;

    // Event interface
public:
    void execute()
    {
        if (cycle() % m_interval == 0)
        {
            m_h5group[to_string(cycle())] = solver().heights();
        }
    }
};

class StoreParticles : public SOSEvent
{
public:
    StoreParticles(const SOSSolver &solver,
                   const OfflatticeMonteCarlo &omc,
                   const uint interval,
                   H5Wrapper::Member &h5group) :
        SOSEvent(solver, "storeparticles"),
        m_omc(omc),
        m_interval(interval),
        m_h5group(h5group.addMember("stored_particles"))
    {

    }

private:
    const OfflatticeMonteCarlo &m_omc;
    const uint m_interval;
    H5Wrapper::Member &m_h5group;

    // Event interface
public:
    void execute()
    {
        if (cycle() % m_interval == 0)
        {
            if (m_omc.nOfflatticeParticles() == 0)
            {
                m_h5group[to_string(cycle())] = mat(3, 0);
            }
            else
            {
                m_h5group[to_string(cycle())] = m_omc.particlePositions()(span::all,
                                                                          span(0, m_omc.nOfflatticeParticles() - 1)).eval();
            }
        }
    }
};
