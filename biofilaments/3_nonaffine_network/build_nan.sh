find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_nan.for 
gfortran -o nan.o *.for *.f90
gfortran -o nan.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_nan.for 
