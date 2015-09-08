#include "firstpassagecontinuum.h"

FirstPassageContinuum::FirstPassageContinuum(SOSSolver &solver,
                                             const double maxdt,
                                             rate_func_type rateLaw) :
    Diffusion(solver, "FirstPassage"),
    OfflatticeMonteCarlo(solver, maxdt),
    m_rateLaw(rateLaw)
{
    insertParticle(0, 0, 2);
    rateLaw(this, 0, 0, 0);
    cout << "....!" << endl;
}

FirstPassageContinuum::~FirstPassageContinuum()
{

}

void FirstPassageContinuum::execute()
{
    //pass
}
