TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../app_defaults.pri)

SOURCES += main.cpp \
    solidonsolidsolver.cpp \
    solidonsolidevents.cpp

QMAKE_CXXFLAGS += -std=c++11

QMAKE_LIBS += -larmadillo

HEADERS += \
    solidonsolidsolver.h \
    solidonsolidevents.h


OTHER_FILES += \
    quasidiffusion.cpp \
    quasidiffusionsystem.cpp \
    movingwall.cpp \
    quasidiffusionevents.cpp \
    quasidiffusion.h \
    quasidiffusionsystem.h \
    movingwall.h \
    quasidiffusionevents.h
