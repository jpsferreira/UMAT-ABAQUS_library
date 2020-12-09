find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_og.for 
gfortran -Wextra  -pedantic -o og *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_og.for 
