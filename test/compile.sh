cd ../src
gfortran ./cfs.f90     -fPIC -larpack -fopenmp -c 
gfortran ./bs.f90      -fPIC -larpack -fopenmp -c 
gfortran ./op.f90      -fPIC -larpack -fopenmp -c 
gfortran ./diag.f90    -fPIC -larpack -fopenmp -c 
gfortran ./diag_re.f90 -fPIC -larpack -fopenmp -c 
gfortran ./ent.f90     -fPIC -larpack -fopenmp -c 
gfortran ./*.o -shared -fPIC -larpack -fopenmp -o libfuzzified.so 
cd ../test