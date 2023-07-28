#!/bin/bash

root=$PWD

# generate libgulpklmc (archiving)
cd ./../Linux_MPI
  ar rcv libgulpklmc.a $( ls *.o | grep -v gulp.o )

cd $root
#ftn -c gulpklmc.F90 -DKLMC -I./../Linux_MPI
#ftn -o interface_test.x gulpklmc.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF

mpicc -c c_main.c
mpif90 -c gulpklmc.F90 -I./../Linux_MPI
mpif90 -c gulpklmc_initmax.F90 -I./../Linux_MPI
mpicc -o klmc.x c_main.o gulpklmc.o gulpklmc_initmax.o -L./../Linux_MPI -lgulpklmc -lifcore -lifport -L./../../Utils/pGFNFF/Src -lpGFNFF -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -lmpi
#mpicc -o klmc.x c_main.o gulpklmc.o gulpklmc_initmax.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -lmpi
