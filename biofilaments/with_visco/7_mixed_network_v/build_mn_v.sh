find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_mn_v.for 
gfortran -Wextra  -pedantic -o mn_v *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_mn_v.for 
  