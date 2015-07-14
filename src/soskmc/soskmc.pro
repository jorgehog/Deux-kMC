include(../../defaults.pri)

TEMPLATE = lib

TARGET = ../../lib/SOSkMC

SOURCES += solidonsolidsolver.cpp \
    solidonsolidreaction.cpp \
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
    Events/diffusion/offlatticemontecarloboundary.cpp \
    Events/diffusion/offlatticemontecarlonoboundary.cpp \
    Events/diffusion/offlatticemontecarlo.cpp

HEADERS += \
    solidonsolidsolver.h \
    solidonsolidreaction.h \
    Events/solidonsolidevent.h \
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
    Events/diffusion/offlatticemontecarloboundary.h \
    Events/diffusion/offlatticemontecarlonoboundary.h \
    Events/diffusion/offlatticemontecarlo.h

LIBS += -lkMC

