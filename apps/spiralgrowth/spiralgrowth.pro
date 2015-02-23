TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../app_defaults.pri)

SOURCES += main.cpp \
    solidonsolidsolver.cpp \
    solidonsolidreaction.cpp \
    pressurewall.cpp \
    equilibriater.cpp \
    eqmu.cpp \
    miscevents.cpp

QMAKE_CXXFLAGS += -std=c++11

QMAKE_LIBS += -larmadillo

HEADERS += \
    solidonsolidsolver.h \
    solidonsolidreaction.h \
    pressurewall.h \
    solidonsolidevent.h \
    equilibriater.h \
    eqmu.h \
    miscevents.h


OTHER_FILES += \
    quasidiffusion.cpp \
    quasidiffusionsystem.cpp \
    movingwall.cpp \
    quasidiffusionevents.cpp \
    quasidiffusion.h \
    quasidiffusionsystem.h \
    movingwall.h \
    quasidiffusionevents.h
