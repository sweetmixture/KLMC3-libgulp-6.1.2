### KLMC3-libgulp-6.1.2
  
Libridised version of GULP 6.1.2  
  
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
