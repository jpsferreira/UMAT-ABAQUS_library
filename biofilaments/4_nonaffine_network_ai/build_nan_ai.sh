find . -type f -path '*src/*' -name '_global.*' -exec cat {} +> umat_nan_ai.f90
find . -type f -path '*src/*'  -name '*.f90' -not -name '_global.f90'  -exec cat {} >> umat_nan_ai.f90 \;
gfortran -o nan_ai.o *.f90
find . -type f -path '*src/*' -name '*.f90' -not -name 'GETOUTDIR.f90' -exec cat {} +> test_in_abaqus/umat_nan_ai.f