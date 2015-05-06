include(../../defaults.pri)

TEMPLATE = lib

TARGET = ../../lib/SOSkMC

SOURCES += solidonsolidsolver.cpp \
    solidonsolidreaction.cpp \
    Events/pressurewall.cpp \
    Events/equilibriater.cpp \
    Events/eqmu.cpp \
    Events/miscevents.cpp

HEADERS += \
    solidonsolidsolver.h \
    solidonsolidreaction.h \
    Events/pressurewall.h \
    Events/solidonsolidevent.h \
    Events/equilibriater.h \
    Events/eqmu.h \
    Events/miscevents.h

LIBS += -lkMC

!equals(PWD, $${OUT_PWD}) {
    QMAKE_POST_LINK += $(COPY_DIR) $$OUT_PWD/../lib $$TOP_PWD
}



