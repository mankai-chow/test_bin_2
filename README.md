This is the Fortran source code in the Julia package [FuzzifiED](https://github.com/mankai-chow/FuzzifiED.jl) intended for exact diagonalisation (ED) calculations on fuzzy sphere. 

If you want to compile this code yourself for the use of FuzzifiED, use the following commands
```bash
cd FuzzifiED_Fortran/src/
FFLAGS=(-O3 -fPIC -fopenmp)
for src in cfs.f90 bs.f90 op.f90 diag.f90 diag_re.f90 ent.f90; do
    gfortran "${FFLAGS[@]}" -c ./${src}
done
gfortran "${FFLAGS[@]}" -shared -o "${libdir}/libfuzzified.${dlext}" ./*.o -larpack
cd super
for src in scfs.f90 sbs.f90 sop.f90; do
    gfortran "${FFLAGS[@]}" -c ./${src}
done
gfortran "${FFLAGS[@]}" -shared -o "${libdir}/libfuzzifino.${dlext}" ./*.o
```
If you use 64-bits version of Arpack, you might also need to include `-fdefault-integer-8` in `FFLAGS`. Before running the code, you need to export `${destdir}` as the destination directory, and `$dlext` is either `so`, `dylib` or `ddl` depending on your system. _E.g_, a Linux user who want to store the shared libraries at `~/FuzzifiED_libs` can set
```bash
export dlext=".so"
export destdir="~/FuzzifiED_libs"
```
When running the Julia code of FuzzifiED, the directory of the compiled dynamic library need to be stored into the environment variables `FuzzifiED.Libpath` and `FuzzifiED.Fuzzifino.Libpathino`. In the example above,
```julia
FuzzifiED.Libpath = "~/FuzzifiED_libs/libfuzzified.so"
FuzzifiED.Fuzzifino.Libpathino = "~/FuzzifiED_libs/libfuzzifino.so"
```

This package is developped by Zheng Zhou (周正) and contributors. If you have any questions, please contact [physics@zhengzhou.page](mailto:physics@zhengzhou.page).
