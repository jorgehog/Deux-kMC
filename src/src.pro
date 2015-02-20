include(../defaults.pri)

TEMPLATE = lib

TARGET = ../lib/kMC

HEADERS = RNG/kmcrng.h \
    RNG/rng.h \
    kmcsolver.h \
    kmcassets.h \
    reaction.h \

SOURCES += \
    RNG/rng.cpp \
    kmcsolver.cpp \
    reaction.cpp \

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

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$TOP_PWD
}



