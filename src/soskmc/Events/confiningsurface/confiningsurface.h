#pragma once

#include "../sosevent.h"
#include "../../heightconnecter.h"

class DissolutionDeposition;

class ConfiningSurface : public SOSEvent, public HeightConnecter
{
public:
    ConfiningSurface(SOSSolver &solver,
                     string type,
                     string unit = "",
                     bool hasOutput = false,
                     bool storeValue = false);

    virtual ~ConfiningSurface();

    virtual bool hasSurface() const = 0;

    virtual double confinementEnergy(const uint x, const uint y) = 0;

    virtual bool acceptDiffusionMove(const double x0, const double y0, const double z0,
                                     const double x1, const double y1, const double z1) const = 0;

    virtual double diffusionDrift(const double x, const double y, const double z) const = 0;

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
