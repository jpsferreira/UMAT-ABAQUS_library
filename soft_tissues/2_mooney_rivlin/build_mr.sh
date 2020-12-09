find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_mr.for 
gfortran -Wextra  -pedantic -o mr *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_mr.for 
