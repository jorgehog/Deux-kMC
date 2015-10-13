TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt

SUBDIRS += ignis \
    HDF5Wrapper \
    intercombinatorzor

HEADERS += lammpswriter/lammpswriter.h \
    libconfig_utils/libconfig_utils.h \
    BADAss/badass.h

OTHER_FILES += utils.h
