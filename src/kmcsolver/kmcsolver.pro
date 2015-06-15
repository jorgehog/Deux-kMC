include(../../defaults.pri)

TEMPLATE = lib

TARGET = ../../lib/kMC

HEADERS = RNG/kmcrng.h \
    RNG/rng.h \
    kmcsolver.h \
    kmcassets.h \
    reaction.h \
    boundary/boundary.cpp \
    boundary/edge.cpp \
    boundary/periodic.cpp

SOURCES += \
    RNG/rng.cpp \
    reaction.cpp \
    kmcsolver.cpp \
    boundary/boundary.cpp \
    boundary/edge.cpp \
    boundary/periodic.cpp

RNG_ZIG
{

DEFINES += KMC_RNG_ZIG

HEADERS += RNG/kmcrngzig.h \
           RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}

QMAKE_PRE_LINK += $(MKDIR) $$PWD/../lib $$shadowed($$PWD)/../lib

