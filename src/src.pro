TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += kmcsolver \
    soskmc

soskmc.depends = kmcsolver


