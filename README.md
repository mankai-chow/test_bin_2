This is the Fortran source code in the Julia package [FuzzifiED](https://github.com/FuzzifiED/FuzzifiED.jl) intended for exact diagonalisation (ED) calculations on fuzzy sphere. 

If you want to compile this code yourself for the use of FuzzifiED, please refer to `build_here.sh`. If you use 32-bits version of Arpack, remove  `-fdefault-integer-8` from `FFLAGS`. If you are not using Linux, change `.so` to `.dylib` or `.ddl`. When running the Julia code of FuzzifiED, the directory of the compiled dynamic library need to be stored into the environment variables `FuzzifiED.Libpath` and `FuzzifiED.Fuzzifino.Libpathino`. 

This package is developed by Zheng Zhou (周正) and contributors. If you have any questions, please contact [fuzzified@zhengzhou.page](mailto:fuzzified@zhengzhou.page).
