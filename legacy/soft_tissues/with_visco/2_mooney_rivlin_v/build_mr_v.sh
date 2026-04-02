find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_mr_v.for 
gfortran -o mr_v.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_mr_v.for 
