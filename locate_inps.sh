cwd=${PWD##*/}
for file in $(find . -name 'cube_*.inp' -type f); do
echo $file
echo $(basename $file)
done