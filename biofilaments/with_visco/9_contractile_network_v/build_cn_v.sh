find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_cn_v.for 
gfortran -o cn_v.o *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_cn_v.for 
  