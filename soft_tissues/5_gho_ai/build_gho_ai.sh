find . -type f -path '*src/*' -name '_global.*' -exec cat {} +> umat_gho_ai.f90
find . -type f -path '*src/*'  -name '*.f90' -not -name '_global.f90'  -exec cat {} >> umat_gho_ai.f90 \;
gfortran -Wextra  -pedantic -o gho_ai *.f90
gfortran -Wextra  -pedantic -o gho_ai *.f90
find . -type f -path '*src/*' -name '_global.*' -exec cat {} +> test_in_abaqus/umat_gho_ai.f
find . -type f -path '*src/*'  -name '*.f90' -not -name 'GETOUTDIR.f90' -not -name '_global.f90'  -exec cat {} >> test_in_abaqus/umat_gho_ai.f \;
