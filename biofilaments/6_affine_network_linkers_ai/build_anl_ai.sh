#ls -rd *src*/*.f90 | LC_ALL=C sort | xargs cat > umat_anl_ai.f90
find . -type f -path '*src/*' -name '_global.*' -exec cat {} +> umat_anl_ai.f90
find . -type f -path '*src/*'  -name '*.f90' -not -name '_global.f90'  -exec cat {} >> umat_anl_ai.f90 \;
gfortran -o anl_ai *.f90
find . -type f -path '*src/*' -name '*.f90' -not -name 'GETOUTDIR.f90' -exec cat {} +> test_in_abaqus/umat_anl_ai.f