#pragma once

#include "sosevent.h"


class EqMu : public SOSEvent
{
public:

    EqMu(const SOSSolver &solver) :
        SOSEvent(solver, "EqMu", "", true, true),
        m_accuNeighbours(0),
        m_accuDissolutionRate(0),
        m_totalTime(0)
    {

    }

    ~EqMu()
    {

    }

    void initialize();

    void execute();

    void reset();

    void restart();

    double dMu() const
    {
        return log(m_dMu);
    }

private:

    double m_dMu;

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

