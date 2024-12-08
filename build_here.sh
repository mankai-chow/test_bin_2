mkdir lib_local
cd src/
FFLAGS=(-O3 -fPIC -fopenmp -fdefault-integer-8)
for src in cfs.f90 bs.f90 op.f90 diag.f90 diag_re.f90 ent.f90; do
    gfortran "${FFLAGS[@]}" -c ./${src}
done
gfortran "${FFLAGS[@]}" -shared -o "../lib_local/libfuzzified.so" ./*.o -larpack
cd super
for src in scfs.f90 sbs.f90 sop.f90 sent.f90; do
    gfortran "${FFLAGS[@]}" -c ./${src}
done
gfortran "${FFLAGS[@]}" -shared -o "../../lib_local/libfuzzifino.so" ./*.o

rm *.mod *.o 
cd ..
rm *.mod *.o 
cd ..