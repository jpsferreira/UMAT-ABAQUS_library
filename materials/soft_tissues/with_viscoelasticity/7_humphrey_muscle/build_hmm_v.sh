find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_hmm_v.for 
gfortran -Wextra  -pedantic -o hmm_v *.for *.f90
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_hmm_v.for 
