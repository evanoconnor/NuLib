#!/bin/bash
module load f90wrap
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/ceclub/gr1dnulib/mesasdk/lib:$(readlink -e .)
