find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_hm_d.for 
gfortran -Wextra  -pedantic -o hm_d *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_hm_d.for 
