#!/bin/bash
# ---------------------------------------------------------------------------
#
#   Author:         Woongkyu Jee / woong.jee.16@ucl.ac.uk
#   Affiliation:    University College London
#   Date:           2023.05.25 - 
#
# ---------------------------------------------------------------------------

root=$PWD
exe="call_gulpmain.x"

# ---------------------------------------------------------------------------
# generate libgulpklmc (archiving)
#
# before run this script, build gulp first at ${GULPROOT}/Src
# e.g., $ mkgulp_intel -c intel -j 4 -m
# -c : compiler type
# -j : parallel make
# -m : using MPI
# ---------------------------------------------------------------------------
cd ./../Linux_MPI
  ar rcv libgulpklmc.a $( ls *.o | grep -v gulp.o )
cd $root

# ---------------------------------------------------------------------------
# wkjee 06.2024
# standard compilation (on Tier-2 HPCs, MMM YOUNG)
#
# tested environement?
# ---------------------------------------------------------------------------

mpicc  -c call_gulpmain.c
mpif90 -c gulpklmc.F90 -I./../Linux_MPI
mpif90 -c gulpklmc_initmax.F90 -I./../Linux_MPI
mpif90 -c gulpklmc_deallocate_all.F90 -I./../Linux_MPI
mpicc  -o ${exe} call_gulpmain.o gulpklmc.o gulpklmc_initmax.o gulpklmc_deallocate_all.o \
          -L./../Linux_MPI -lgulpklmc -lifcore -lifport \
          -L./../../Utils/pGFNFF/Src -lpGFNFF \
          -L${MKLROOT}/lib/intel64 \
          -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl -lmpi
rm *.o

# ---------------------------------------------------------------------------
# * * * Script Outcome * * *
# After running this script, gulplib 6.1.2 will be generated in 'Linux_MPI' (libgulpklmc.a)
# ---------------------------------------------------------------------------
