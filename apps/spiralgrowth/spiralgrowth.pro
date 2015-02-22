TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../app_defaults.pri)

SOURCES += main.cpp \
    spiralgrowthsolver.cpp

QMAKE_CXXFLAGS += -std=c++11

QMAKE_LIBS += -larmadillo

HEADERS += \
    spiralgrowthsolver.h


OTHER_FILES += \
    quasidiffusion.cpp \
    quasidiffusionsystem.cpp \
    movingwall.cpp \
    quasidiffusion.h \
    quasidiffusionsystem.h \
    movingwall.h
