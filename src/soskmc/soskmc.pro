include(../../defaults.pri)

TEMPLATE = lib

TARGET = ../../lib/SOSkMC

SOURCES += solidonsolidsolver.cpp \
    solidonsolidreaction.cpp \
    Events/equilibriater.cpp \
    Events/eqmu.cpp \
    Events/miscevents.cpp \
    Events/cavitydiffusion.cpp \
    Events/confiningsurface/confiningsurface.cpp \
    Events/confiningsurface/fixedrdlsurface.cpp \
    Events/confiningsurface/rdlsurface.cpp \
    Events/confiningsurface/fixedsurface.cpp

HEADERS += \
    solidonsolidsolver.h \
    solidonsolidreaction.h \
    Events/solidonsolidevent.h \
    Events/equilibriater.h \
    Events/eqmu.h \
    Events/miscevents.h \
    Events/cavitydiffusion.h \
    Events/confiningsurface/confiningsurface.h \
    Events/confiningsurface/fixedrdlsurface.h \
    Events/confiningsurface/rdlsurface.h \
    Events/confiningsurface/fixedsurface.h

LIBS += -lkMC

