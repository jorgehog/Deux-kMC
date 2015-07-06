#ifndef CONFININGSURFACE_H
#define CONFININGSURFACE_H

#include "../solidonsolidevent.h"

class DiffusionDeposition;

class ConfiningSurface : public SolidOnSolidEvent
{
public:
    ConfiningSurface(SolidOnSolidSolver &solver,
                     string type,
                     string unit = "",
                     bool hasOutput = false,
                     bool storeValue = false);

    virtual ~ConfiningSurface();

    virtual void setupInitialConditions() = 0;

    virtual void registerHeightChange(const uint x,
                                      const uint y,
                                      std::vector<DiffusionDeposition *> affectedReactions,
                                      const uint n) = 0;

    virtual double confinementEnergy(const uint x, const uint y) = 0;

    const double &height() const
    {
        return m_height;
    }

    void setHeight(const double height)
    {
        m_height = height;
    }

private:

    double m_height;
};

#endif // CONFININGSURFACE_H
