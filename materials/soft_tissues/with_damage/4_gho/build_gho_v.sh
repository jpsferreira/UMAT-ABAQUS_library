find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_gho_d.for 
gfortran -Wextra  -pedantic -o gho_d *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_gho_d.for 
