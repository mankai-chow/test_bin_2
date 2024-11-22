cd ../src
gfortran -fPIC -c ./cfs.f90      
gfortran -fPIC -c ./bs.f90       
gfortran -fPIC -c ./op.f90       
gfortran -fPIC -c ./diag.f90     
gfortran -fPIC -c ./diag_re.f90  
gfortran -fPIC -c ./ent.f90      
gfortran ./*.o -shared -larpack -fopenmp -o libfuzzified.so 
cd ../test