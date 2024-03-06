#!/usr/bin/env python 

""" 
    Script that takes a NuLib table, checks
    if the source code is present. If so, 
    it extracts the source tarball from 
    the hdf5 file and places it into a 
    separate directory.

    Requires h5py and numpy and their dependencies.
"""

import os
import sys
import fnmatch
import glob
import h5py
import numpy as np


if len(sys.argv) < 2:
    print "Usage: nulib_extract_source_in_table.py <NuLib Table Name>"
    sys.exit()

nulib_table_name = sys.argv[1]

tarfile = "nulib_src.tar.gz"

# creating output directory
print "Creating output directory saved_nulib"
os.system("mkdir saved_nulib")

h5file = h5py.File(nulib_table_name,"r")

try:
    indata = h5file['NuLib Source'][()]
except:
    print "Sorry, no source availabe in this NuLib file."
    h5file.close()
    sys.exit()

h5file.close()

outfile = open("saved_nulib/"+tarfile,"wb")
outfile.write(indata)
outfile.close()


