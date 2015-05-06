include(../../defaults.pri)

TEMPLATE = lib

TARGET = ../../lib/kMC

HEADERS = RNG/kmcrng.h \
    RNG/rng.h \
    kmcsolver.h \
    kmcassets.h \
    reaction.h \

SOURCES += \
    RNG/rng.cpp \
    reaction.cpp \
    kmcsolver.cpp

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


