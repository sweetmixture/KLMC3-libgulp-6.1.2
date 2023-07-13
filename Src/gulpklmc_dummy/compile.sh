#!/bin/bash

root=$PWD

# generate libgulpklmc (archiving)
cd ./../Linux_MPI
  ar rcv libgulpklmc.a $( ls *.o | grep -v gulp.o )

cd $root
#ftn -c gulpklmc.F90 -DKLMC -I./../Linux_MPI
#ftn -o interface_test.x gulpklmc.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF

mpicc -c c_main.c
ftn -c gulpklmc.F90 -I./../Linux_MPI
ftn -o klmc.x c_main.o gulpklmc.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF
