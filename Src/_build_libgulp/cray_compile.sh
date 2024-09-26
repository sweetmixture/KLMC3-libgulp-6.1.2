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
# e.g. 1 using gnu, $ mkgulp_cray_gnu -c cray -j 4 -m
# e.g. 2 using cce, $ mkgulp_cray_cce -c cray -j 4 -m
# -c : compiler type
# -j : parallel make
# -m : using MPI
# ---------------------------------------------------------------------------
cd ./../Linux_MPI
  ar rcv libgulpklmc.a $( ls *.o | grep -v gulp.o )
cd $root

# ---------------------------------------------------------------------------
# wkjee 08.2023
# -fallow-argument-mismatch : ignoring type mismathcing on using gnu compiler (cray)
# 
# standard compilation (on ARCHER 2, using 'cce' and 'gnu' stable)
# ---------------------------------------------------------------------------
cc -c call_gulpmain.c
ftn -fallow-argument-mismatch -c gulpklmc.F90 -I./../Linux_MPI
ftn -fallow-argument-mismatch -c gulpklmc_initmax.F90 -I./../Linux_MPI
ftn -fallow-argument-mismatch -c gulpklmc_deallocate_all.F90 -I./../Linux_MPI
ftn -o ${exe} call_gulpmain.o gulpklmc.o gulpklmc_initmax.o gulpklmc_deallocate_all.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF
rm *.o

# ---------------------------------------------------------------------------
# wkjee 06.2024 : for development only (using SCORE-P debugger)
#                 to fully utilised SCORE-P functionality,
#                 please compile 'gulp' with 'mkgulp_scorep' in ${GULPROOT}/Src
#                 e.g., $ ./mkgulp_scorep -c cray -j 4 -m
# 
# 06.2024: memory leak fix (force to dealloac all) -> see 'klmcgulp_lib.log'
w
#          to use scorep compiler, set relevant modules
#          e.g., on ARCHER 2
#          module unload perftools-base/22.12.0
#          module load other-software
#          module load PrgEnv-gnu
#          module load scalasca/2.6.1-gcc11
#
# ---------------------------------------------------------------------------
#scorep-cc -c call_gulpmain.c
#scorep-ftn -fallow-argument-mismatch -c gulpklmc.F90 -I./../Linux_MPI
#scorep-ftn -fallow-argument-mismatch -c gulpklmc_initmax.F90 -I./../Linux_MPI
#scorep-ftn -fallow-argument-mismatch -c gulpklmc_deallocate_all.F90 -I./../Linux_MPI
#scorep-ftn -o ${exe} call_gulpmain.o gulpklmc.o gulpklmc_initmax.o gulpklmc_deallocate_all.o -L./../Linux_MPI -lgulpklmc -L./../../Utils/pGFNFF/Src -lpGFNFF
#rm *.o


# ---------------------------------------------------------------------------
# * * * Script Outcome * * *
# After running this script, gulplib 6.1.2 will be generated in 'Linux_MPI' (libgulpklmc.a)
# ---------------------------------------------------------------------------
