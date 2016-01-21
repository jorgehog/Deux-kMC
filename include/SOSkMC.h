#pragma once

#include "kMC.h"

static_assert(sizeof(arma::sword) == 4, "64 bit armadillo is not supported.");

#include "../src/soskmc/sossolver.h"

#include "../src/soskmc/Events/confiningsurface/noconfinement.h"
#include "../src/soskmc/Events/confiningsurface/constantconfinement.h"
#include "../src/soskmc/Events/confiningsurface/fixedsurface.h"
#include "../src/soskmc/Events/confiningsurface/rdlsurface.h"
#include "../src/soskmc/Events/confiningsurface/fixedrdlsurface.h"

#include "../src/soskmc/Events/confiningsurface/rdlpotential.h"

#include "../src/soskmc/Events/diffusion/constantconcentration.h"
#include "../src/soskmc/Events/diffusion/latticediffusion.h"
#include "../src/soskmc/Events/diffusion/fixedpointtimestepping.h"
#include "../src/soskmc/Events/diffusion/astarfirstpassage.h"
#include "../src/soskmc/Events/diffusion/radialfirstpassage.h"
#include "../src/soskmc/Events/diffusion/multiscale.h"
#include "../src/soskmc/Events/diffusion/confinedconstantconcentration.h"

#include "../src/soskmc/Events/eqmu.h"
#include "../src/soskmc/Events/equilibriater.h"
#include "../src/soskmc/Events/miscevents.h"

#include "../src/soskmc/Events/dumpsystem.h"

#include "../src/soskmc/averageheightboundary.h"

#include "../src/soskmc/sosdiffusionreaction.h"
#include "../src/soskmc/dissolutiondeposition.h"

#include "../src/soskmc/particlenumberconservator.h"
#include "../src/soskmc/localpotential.h"
#include "../src/soskmc/localcachedpotential.h"
