find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_an_v.for 
gfortran -o an_v.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_an_v.for 
