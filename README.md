### General Utility Lattice Program (GULP) Static Library for Knowledge Led Master Code 3 (KLMC3)
---
Original GULP source code is from: https://gulp.curtin.edu.au.

This modification has been made to use ```gulpmain()``` soubroutine, for instance, which could be interfaced by using a wrapper (currently used in KLMC3: https://github.com/sweetmixture/KLMC3).


  * Examples of building GULP static library

    As a common step, clone this git-repo on your Linux system
    ```
    git clone https://github.com/sweetmixture/KLMC3-libgulp-6.1.2.git
    ```
    and move to ```/KLMC3-libgulp-6.1.2.git/Src/```.

    Build of the library has been tested on CRAY (GNU) and INTEL systems.

    ---
    #### [A] On CRAY system
    for instance, if you have an account on ARCHER 2 (https://www.archer2.ac.uk):

    STEP1. First make sure to use following modules,
    ```
    $ module restore
    $ module load PrgEnv-gnu
    ```
    module ```PrgEnv-gnu``` could be replaced to any other equivalents, e.g., ```PrgEnv-gnu-amd```.

    STEP2. Creating GULP *.mod / *.o files,
    ```
    $ cd /path/to/KLMC3-libgulp-6.1.2/Src
    $ ./mkgulp_cray_gnu -c cray -j 4 -m
    ```
    where  
      -c cray : using "cray" compiler  
      -j 4    : using 4 cpu cores for this makefile run  
      -m      : using MPI  

    STEP3. Creating static library  
    ```
    cd _build_libgulp
    bash cray_compile.sh
    ```
    This will generate ```libgulpklmc.a``` at: ```/path/to/KLMC3-libgulp-6.1.2/Src/Linux_MPI```,
    which will be automatically picked up by KLMC3 during its compilation.

    ---
    #### [B] On Intel system
    for instance, if you have an account on UCL MMM YOUNG (https://www.rc.ucl.ac.uk/docs/Clusters/Young):

    STEP1. First make sure to use following modules,
    ```
    $ module load gcc-libs/4.9.2
    $ module load compilers/intel/2019/update4
    $ module load mpi/intel/2019/update4/intel
    ```
    It could be possible to use different versions of INTEL compilers and INTEL-MPI libraries, however, other verions have not checked.  

    * The rest of steps (2 and 3) are as in [A] with slightly different commands.  
   
    STEP2. Creating GULP *.mod / *.o files,
    ```
    $ cd /path/to/KLMC3-libgulp-6.1.2/Src
    $ ./mkgulp_intel -c intel -j 4 -m
    ```
    where  
      -c intel : using "intel" compiler  
      -j 4    : using 4 cpu cores for this makefile run  
      -m      : using MPI  

    STEP3. Creating static library  
    ```
    cd _build_libgulp
    bash cray_compile.sh
    ```
    This will generate ```libgulpklmc.a``` at: ```/path/to/KLMC3-libgulp-6.1.2/Src/Linux_MPI```,
    which will be automatically picked up by KLMC3 during its compilation.
    
    
  
