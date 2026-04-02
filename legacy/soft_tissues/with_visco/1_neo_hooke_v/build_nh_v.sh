find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_nh_v.for 
gfortran -o nh_v.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_nh_v.for 
