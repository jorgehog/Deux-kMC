#include "firstpassagecontinuum.h"

FirstPassageContinuum::FirstPassageContinuum(SOSSolver &solver,
                                             const double maxdt,
                                             rate_func_type rateLaw) :
    Diffusion(solver, "FirstPassage"),
    OfflatticeMonteCarlo(solver, maxdt),
    m_rateLaw(rateLaw)
{

}

FirstPassageContinuum::~FirstPassageContinuum()
{

}

void FirstPassageContinuum::execute()
{
    //pass
}
