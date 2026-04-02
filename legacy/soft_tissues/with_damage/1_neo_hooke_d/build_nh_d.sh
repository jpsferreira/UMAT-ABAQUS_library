find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_nh_d.for 
gfortran -o nh_d.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_nh_d.for 
