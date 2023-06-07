find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_mn.for 
gfortran -o mn.o *.for *.f90
gfortran -o mn.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_mn.for 
  