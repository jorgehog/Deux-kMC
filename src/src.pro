include(../defaults.pri)

TEMPLATE = lib

TARGET = ../lib/kMC

HEADERS = RNG/kMCRNG.h

SOURCES +=

RNG_ZIG {

HEADERS += RNG/zigrandom.h \
           RNG/zignor.h

SOURCES += RNG/zigrandom.cpp \
           RNG/zignor.cpp

}

QMAKE_PRE_LINK += $(MKDIR) $$PWD/../lib $$shadowed($$PWD)/../lib

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$TOP_PWD
}



