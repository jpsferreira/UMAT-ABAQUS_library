find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_gho.for 
gfortran -Wextra  -pedantic -o gho *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_gho.for 
