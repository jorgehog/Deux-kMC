#pragma once

#include "kmcrng.h"

#include "zignor.h"
#include "zigrandom.h"


namespace kMC
{

class KMCRNGZIG : public KMCRNG<int>
{
public:
    inline double normal();

    inline double uniform();

private:
    inline void onInitialize();

};

void KMCRNGZIG::onInitialize()
{
    int inseed = static_cast<int>(initialSeed());
    int cseed = 100;
    int seed2 = inseed * 3;
    RanSetSeed_MWC8222(&seed2, cseed);
    RanNormalSetSeedZig32(&inseed, 5);
}

double KMCRNGZIG::normal()
{
    return DRanNormalZig32();
}

double KMCRNGZIG::uniform()
{
    return DRan_MWC8222();
}

}
