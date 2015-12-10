#include "concentrationprofile.h"

#include "../../sossolver.h"


ConcentrationProfile::ConcentrationProfile(SOSSolver &solver, funcType profileFunction) :
    Diffusion(solver, "ConcProfile"),
    m_profileFunction(profileFunction)
{

}

ConcentrationProfile::~ConcentrationProfile()
{

}

double ConcentrationProfile::depositionRate(const uint x, const uint y) const
{

    //since profile function gives dep rate in unit concentration, we have to check for this.
    if (solver().concentrationIsZero())
    {
        return 0;
    }

    return m_profileFunction(x, y);
}

uint ConcentrationProfile::dissolutionPaths(const uint x, const uint y) const
{
    (void) x;
    (void) y;

    return 1u;
}

void ConcentrationProfile::notifyObserver(const Subjects &subject)
{
    (void) subject;
}

void ConcentrationProfile::executeDiffusionReaction(SOSDiffusionReaction *reaction, const int x, const int y, const int z)
{
    (void) reaction;
    (void) x;
    (void) y;
    (void) z;
}

void ConcentrationProfile::executeConcentrationBoundaryReaction(ConcentrationBoundaryReaction *reaction)
{
    (void) reaction;
}

bool ConcentrationProfile::isBlockedPosition(const uint x, const uint y, const int z) const
{
    (void) x;
    (void) y;
    (void) z;

    return false;
}

double ConcentrationProfile::concentration() const
{
    return solver().concentration();
}

bool ConcentrationProfile::hasDiscreteParticles() const
{
    return false;
}

