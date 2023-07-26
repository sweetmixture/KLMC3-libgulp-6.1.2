#!/bin/bash
ifort -c fortsub.f90
icc -c c_main.c
icc -o a.out c_main.o fortsub.o -lifcore -limf
