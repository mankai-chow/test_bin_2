This is the Fortran source code in the Julia package [FuzzifiED](https://github.com/mankai-chow/FuzzifiED.jl) intended for exact diagonalisation (ED) calculations on fuzzy sphere. 

If you want to compile this code yourself for the use of FuzzifiED, use the following commands
```
gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./cfs.f90
gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./bs.f90
gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./op.f90
gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./diag.f90
gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./diag_re.f90
gfortran -fPIC -larpack -fopenmp -L ${libdir} -c ./ent.f90
gfortran -fPIC -shared -larpack -fopenmp -L ${libdir} -o ${destdir}/libfuzzified.$dlext ./*.o
```
where `${libdir}` is the path where the Arpack library is stored, `${destdir}` is the destination directory, and `$dlext` is either `so`, `dylib` or `ddl` depending on your system, and store the directory of the compiled dynamic library into the environment variable `FuzzifiED.Libpath`.

This package is developped by Zheng Zhou (周正) at Perimeter Institute and collaborators. If you have any questions, please contact at zzhou@pitp.ca.
