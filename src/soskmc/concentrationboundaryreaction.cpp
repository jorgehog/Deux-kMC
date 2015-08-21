#include "concentrationboundaryreaction.h"
#include "sossolver.h"
#include "Events/diffusion/diffusion.h"

concentrationBoundaryReaction::concentrationBoundaryReaction(const uint dim, const uint orientation, SOSSolver &solver) :
    Reaction(),
    m_dim(dim),
    m_orientation(orientation),
    m_solver(solver)
{

}

concentrationBoundaryReaction::~concentrationBoundaryReaction()
{

}

bool concentrationBoundaryReaction::isAllowed() const
{
    return !solver().diffusionEvent().isBlockedPosition(x(), y(), z());
}

void concentrationBoundaryReaction::executeAndUpdate()
{
    //choose z randomly across boundary. Subtract Nb and find correct z.

    solver().diffusionEvent().insertDiffusingParticle(x(), y(), z());
}

double concentrationBoundaryReaction::rateExpression()
{
    //find dh and nb.

    return (dh - Nb)/(1-exp(solver().gamma() - 2*solver().alpha()));
}
