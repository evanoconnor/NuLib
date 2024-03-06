#!/usr/bin/env python 

""" 
    Script that takes the NuLib source code and
    parameter file, creates a tarball, and includes
    that tarball in the NuLib table.
    Requires h5py and numpy and their dependencies.
"""

import os
import sys
import fnmatch
import glob
import h5py
import numpy as np


if len(sys.argv) < 2:
    print "Usage: nulib_include_source_in_table.py <NuLib Table Name>"
    sys.exit()

nulib_table_name = sys.argv[1]

# make file list
filelist = []
filelist.append("Makefile")
filelist.append("make.inc")
filelist.append("parameters")
filelist.append("README")

matches = []
for root, dirnames, filenames in os.walk('src'):
    for filename in fnmatch.filter(filenames, '*.F90') + \
        fnmatch.filter(filenames, '*.f') + fnmatch.filter(filenames, 'Makefile') +\
        fnmatch.filter(filenames, 'NuLib_README') + fnmatch.filter(filenames, '*.inc'):
        matches.append(os.path.join(root, filename))

filelist = filelist + matches

tarfile = "nulib_src.tar.gz"
tarstring = "tar -czvf " + tarfile + " "
for xfile in filelist:
    tarstring = tarstring + xfile + " "

print tarstring
os.system(tarstring)

infile = open(tarfile,"rb")
data = infile.read()
infile.close()

wrapdata = np.void(data)

h5file = h5py.File(nulib_table_name,"r+")

dset = h5file.create_dataset("NuLib Source", data=wrapdata)

h5file.close()
