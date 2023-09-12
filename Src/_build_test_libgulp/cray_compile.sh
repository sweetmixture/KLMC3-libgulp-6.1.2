#!/bin/bash

root=$PWD

# generate libgulpklmc (archiving)
cd ./../Linux_MPI
  ar rcv libgulpklmc.a $( ls *.o | grep -v gulp.o )
cd $root

#
# wkjee 08/23
# -fallow-argument-mismatch : ignoring type mismathcing on gnu compiler (cray)
#

mpicc -c c_main.c
ftn -fallow-argument-mismatch -c gulpklmc.F90 -I./../Linux_MPI
ftn -fallow-argument-mismatch -c gulpklmc_initmax.F90 -I./../Linux_MPI
ftn -o klmc.x c_main.o gulpklmc.o gulpklmc_initmax.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF
rm *.o








# below command is not working: why?
#mpicc -o klmc.x c_main.o gulpklmc.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF
