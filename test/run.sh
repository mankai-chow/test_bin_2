gfortran test_diag.f90 -L ../src -I ../src -lfuzzified -o test_diag.out && ./test_diag.out 
gfortran test_diag_re.f90 -L ../src -I ../src -lfuzzified -o test_diag_re.out && ./test_diag_re.out 