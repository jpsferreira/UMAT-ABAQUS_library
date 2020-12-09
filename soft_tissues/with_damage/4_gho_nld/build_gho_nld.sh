find . -type f -path '*src/*' -name '*.for' -exec cat {} +> umat_gho_nld.for 
find . -type f -path '*src/*' -name '*.for' -not -name 'GETOUTDIR.for' -exec cat {} +> test_in_abaqus/umat_gho_nld.for 
