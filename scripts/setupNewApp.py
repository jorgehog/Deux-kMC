import sys, os, re
from os.path import join

if len(sys.argv) == 1:
    raise RuntimeError("Usage: python %s [app name]" % sys.argv[0])

appName = sys.argv[1]
appNameC = appName[0].upper() + appName[1:]

cwd = os.getcwd()
appDir = join(cwd, "..", "apps")
absAppPath = join(appDir, appName)
infilesPath = join(absAppPath, "infiles")

if not os.path.exists(appDir):
    raise RuntimeError("App dir not found: %s does not exist" % appDir)

if os.path.exists(absAppPath):
    raise NameError("App %s already exists." % appName)

os.mkdir(absAppPath)
os.mkdir(infilesPath)

configFileAbsPath = join(infilesPath, "%s.cfg" % appName)
appProFileAbsPath = join(absAppPath, "%s.pro" % appName)
appMainFileAbsPath = join(absAppPath, "%smain.cpp" % appName)


with open(appProFileAbsPath, 'w') as appProFile:

    with open(join(cwd, "app_defaults", "default.pro.bones"), 'r') as defaultFile:

        appProFile.write(defaultFile.read().replace("__name__", appName))

with open(appMainFileAbsPath, 'w') as appMainFile:

    with open(join(cwd, "defaults", "defaultmain.cpp.bones"), 'r') as defaultFile:

        appMainFile.write(defaultFile.read().replace("__name__", appName).replace("__cname__", appNameC))

with open(configFileAbsPath, 'w') as appConfigFile:

    with open(join(cwd, "defaults", "defaultconfig.cfg.bones"), 'r') as defaultFile:

        appConfigFile.write(defaultFile.read().replace("__name__", appName))


appsProFileAbspath = join(cwd, "apps.pro")

with open(appsProFileAbspath, 'r') as appsProFile:

    newAppsProFile = re.sub("(SUBDIRS\s*\+\=+s*.+\n)(\s+)(\w+)\n",
                            "\g<1>\g<2>\g<3> \\\n\g<2>%s\n" % appName,
                            appsProFile.read(), flags=re.DOTALL)

with open(appsProFileAbspath, 'w') as appsProFile:

    appsProFile.write(newAppsProFile)
