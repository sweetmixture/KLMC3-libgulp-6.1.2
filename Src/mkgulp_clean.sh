cd _build_libgulp
bash clean_testrun.sh
rm -rf *.x *.a
cd ..

cd _build_libgulp_mpierr_07.2024
bash clean_testrun.sh
rm -rf *.x *.a
cd ..

rm -rf Linux_MPI
rm -rf ./../Utils/pGFNFF/Src/*.o
rm -rf ./../Utils/pGFNFF/Src/*.mod
rm -rf ./../Utils/pGFNFF/Src/*.a
