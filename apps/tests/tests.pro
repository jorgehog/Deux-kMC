TARGET = kmctests

include(../app_defaults.pri)

SOURCES += main.cpp

LIBS += -lgtest -lgtest_main -lpthread

OTHER_FILES += \
    infiles/tests.cfg

HEADERS += \
    testdiffusion.h \
    kmctestfixture.h

DEFINES += BADASSNOTHROW
