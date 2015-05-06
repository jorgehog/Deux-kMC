TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += utils src apps

src.depends = utils
apps.depends = src utils

OTHER_FILES += infiles/config.cfg \
               .gitignore \
               include/kMC.h \
               include/SOSkMC.h

include(deployment.pri)
qtcAddDeployment()
