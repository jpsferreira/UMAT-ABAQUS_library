find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_anl_v.for 
gfortran -Wextra  -pedantic -o anl_v *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_anl_v.for 
  