#pragma once

#include <sys/types.h>

namespace kMC
{

class Reaction
{
public:

    Reaction();

    virtual ~Reaction() {}

    virtual bool isAllowed() const = 0;

    virtual void executeAndUpdate() = 0;

    virtual double rateExpression() = 0;

    void calculateRate()
    {
        m_rate = rateExpression();
    }

    const double &rate() const
    {
        return m_rate;
    }

private:

    double m_rate;

};

}

