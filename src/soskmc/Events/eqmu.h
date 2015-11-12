#pragma once

#include "sosevent.h"


class EqGamma : public SOSEvent
{
public:

    EqGamma(const SOSSolver &solver) :
        SOSEvent(solver, "EqGamma", "", true, true),
        m_accuNeighbours(0),
        m_accuDissolutionRate(0),
        m_totalTime(0)
    {

    }

    ~EqGamma()
    {

    }

    void initialize();

    void execute();

    void reset();

    void restart();

    double dGamma() const
    {
        return log(m_dGamma);
    }

private:

    double m_dGamma;

    double m_accuNeighbours;
    double m_accuDissolutionRate;

    double m_totalTime;

    void update();

    void resetCounters()
    {
        m_totalTime = 0;

        m_accuNeighbours = 0;
        m_accuDissolutionRate = 0;
    }


};

