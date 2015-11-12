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

    //the reaction can be affected even if it is not
    //directly involved in the executed reaction
    virtual void affectedUpdateRule()
    {
        calculateRate();
    }

    void calculateRate()
    {
        m_rate = rateExpression();
    }

    const double &rate() const
    {
        return m_rate;
    }

    void changeRate(const double newRate)
    {
        m_rate = newRate;
    }

private:

    double m_rate;

};

}

