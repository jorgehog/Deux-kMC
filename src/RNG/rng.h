#pragma once

#include "kmcrng.h"

#ifdef KMC_RNG_ZIG
#include "kmcrngzig.h"
#endif

namespace kMC
{

#ifdef KMC_RNG_ZIG

static KMCRNGZIG rng;

#endif

}
