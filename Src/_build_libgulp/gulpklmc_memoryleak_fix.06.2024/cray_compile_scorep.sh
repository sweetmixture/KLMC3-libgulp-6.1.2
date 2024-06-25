#!/bin/bash

root=$PWD
exe="call_gulpmain.x"
#
# generate libgulpklmc (archiving)
#

cd ./../Linux_MPI
  ar rcv libgulpklmc.a $( ls *.o | grep -v gulp.o )
cd $root
#
# After running this script, gulplib 6.1.2 will be generated in 'Linux_MPI' (libgulpklmc.a)
#

#
# wkjee 08/23
# -fallow-argument-mismatch : ignoring type mismathcing on using gnu compiler (cray)
#

#
# memory fix 06.2024
#
scorep-cc -c c_main.c
scorep-ftn -fallow-argument-mismatch -c gulpklmc.F90 -I./../Linux_MPI
scorep-ftn -fallow-argument-mismatch -c gulpklmc_initmax.F90 -I./../Linux_MPI
scorep-ftn -fallow-argument-mismatch -c gulpklmc_deallocate.F90 -I./../Linux_MPI
scorep-ftn -fallow-argument-mismatch -c gulpklmc_deallocate_all.F90 -I./../Linux_MPI
scorep-ftn -o ${exe} c_main.o gulpklmc.o gulpklmc_initmax.o gulpklmc_deallocate.o gulpklmc_deallocate_all.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF
rm *.o
