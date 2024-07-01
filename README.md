### KLMC3-libgulp-6.1.2
  
GULP-6.1.2 libridised version:

[1] Aiming to use ```gulpmain()``` soubroutine for repeated calls.

  * To use this functionality, first you must compile the GULP source files, which can be achieved by following command line instructions:

    ! Using CRAY system, e.g., on ARCHER2

    [a] Loading following modules,
    ```
    $ module restore
    $ module load PrgEnv-gnu
    ```
    [b] Creating GULP *.mod / *.o files,
    ```
    $ cd /path/to/KLMC3-libgulp-6.1.2/Src
    $ ./mkgulp -c cray -j 4 -m
    ```
    where
      -c cray : using "cray" compiler
      -j 4    : using 4 cpu cores for this makefile run
      -m      : using MPI

    This

    ... KEEP EDITING

    
  * Each gulpmain() call is intended to conduct a standard GULP run:
    (by 06.2024) For the current implementation, an example wrapper for ```gulpmain()``` subroutine is given at: ```/root/Src/_build_libgulp```
    including following source files
      call_gulpmain.c             :
      gulpklmc.c                  :
      gulpklmc_initmax.c          :
      gulpklmc_deallocate_all.c   :


* ------------------------------------------
* OLD DESCRIPTION - 05.2024
  
Purpose: call gulpmain() multiple times for it's use on the KLMC3 task-farm interface.  
  
* Build (Compilation) : environment ARCHER 2 using cray-GNU

1. In ```/root/Src/```, type command,  
```
  $ ./mkgulp -c cray -j 4 -m
```

2. Once step 1 is finished, type commands,
```
  $ cd _build_libgulp
  $ bash cray_compile.sh
```
which will generate a static library, libgulp.a, in the path ```/root/Src/Linux_MPI/```.  
