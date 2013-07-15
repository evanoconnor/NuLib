F90=gfortran
F90FLAGS= -O3 -ffree-line-length-none #-fopenmp

#if you want to take advantage of openmp use this flag
#F90FLAGS= -O3 -openmp

MODINC="-I ./"
#OMP_NUM_THREADS=12
NUCLEI_HEMPEL=1

#You must have a HDF5 version installed with
#the _same_ compiler, this often means compiling
#it by yourself.
#See README file for compiling HDF5

HDF5DIR=/projects/ceclub/gr1dnulib/mesasdk
HDF5INCS=-I$(HDF5DIR)/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5 -lhdf5_fortran -lhdf5 -lz
