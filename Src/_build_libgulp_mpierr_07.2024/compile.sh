ftn -c write_number.F90
cc -c main.c
cc -o main main.o write_number.o -lgfortran
