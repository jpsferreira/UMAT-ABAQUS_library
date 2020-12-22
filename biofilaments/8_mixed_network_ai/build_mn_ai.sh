find -s . -type f -path '*src/*' -name '*.f90' -exec cat {} +> umat_mn_ai.f90
gfortran -Wextra  -pedantic -o mn_ai *.f90
find . -type f -path '*src/*' -name '*.f90' -not -name 'GETOUTDIR.f90' -exec cat {} +> test_in_abaqus/umat_mn_ai.f