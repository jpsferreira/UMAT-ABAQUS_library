find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_anl.for 
gfortran -o anl.o *.for *.f90
gfortran -o anl.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_anl.for 
  