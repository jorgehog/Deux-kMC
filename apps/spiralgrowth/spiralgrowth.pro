TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(../app_defaults.pri)

SOURCES += main.cpp

QMAKE_CXXFLAGS += -std=c++11

QMAKE_LIBS += -larmadillo
