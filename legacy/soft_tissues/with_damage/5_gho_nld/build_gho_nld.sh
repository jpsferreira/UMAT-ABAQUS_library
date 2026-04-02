find . -type f -path '*src/*' -name '*.for' -not -name 'usdfld.for' -exec cat {} +> umat_gho_nld.for
gfortran -o gho_nld.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_gho_nld.for 
