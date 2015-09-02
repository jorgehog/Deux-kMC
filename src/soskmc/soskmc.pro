include(../../defaults.pri)

TEMPLATE = lib

TARGET = ../../lib/SOSkMC

SOURCES += \
    Events/equilibriater.cpp \
    Events/eqmu.cpp \
    Events/miscevents.cpp \
    Events/confiningsurface/confiningsurface.cpp \
    Events/confiningsurface/fixedrdlsurface.cpp \
    Events/confiningsurface/rdlsurface.cpp \
    Events/confiningsurface/fixedsurface.cpp \
    Events/diffusion/diffusion.cpp \
    Events/diffusion/constantconcentration.cpp \
    Events/confiningsurface/noconfinement.cpp \
    Events/diffusion/offlatticemontecarlo.cpp \
    sossolver.cpp \
    sosreaction.cpp \
    sosdiffusionreaction.cpp \
    dissolutiondeposition.cpp \
    Events/diffusion/latticediffusion.cpp \
    concentrationboundaryreaction.cpp \
    Events/dumpsystem.cpp \
    Events/diffusion/fixedpointtimestepping.cpp \
    Events/diffusion/multiscale.cpp \
    Events/diffusion/firstpassagecontinuum.cpp

HEADERS += \
    Events/equilibriater.h \
    Events/eqmu.h \
    Events/miscevents.h \
    Events/confiningsurface/confiningsurface.h \
    Events/confiningsurface/fixedrdlsurface.h \
    Events/confiningsurface/rdlsurface.h \
    Events/confiningsurface/fixedsurface.h \
    Events/diffusion/diffusion.h \
    Events/diffusion/constantconcentration.h \
    Events/confiningsurface/noconfinement.h \
    Events/diffusion/offlatticemontecarlo.h \
    sossolver.h \
    sosreaction.h \
    sosdiffusionreaction.h \
    dissolutiondeposition.h \
    Events/diffusion/latticediffusion.h \
    concentrationboundaryreaction.h \
    Events/dumpsystem.h \
    Events/sosevent.h \
    Events/diffusion/fixedpointtimestepping.h \
    Events/diffusion/multiscale.h \
    Events/diffusion/firstpassagecontinuum.h

LIBS += -lkMC

