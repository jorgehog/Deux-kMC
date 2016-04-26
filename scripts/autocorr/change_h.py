import re
import sys
import glob
import os
import shutil

path = sys.argv[1]
orig = sys.argv[2]
dest = sys.argv[3]

for f_orig in glob.glob(path + "/acf_h*"):
    
    f_dest = re.sub("h%s\_" % orig, "h%s_" % dest, f_orig)
    
    try:
        shutil.copy(f_orig, f_dest)
    except:
        pass

    print "yo, ", f_dest


