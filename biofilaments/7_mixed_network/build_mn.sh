find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_mn.for 
gfortran -Wextra  -pedantic -o mn *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_mn.for 
  