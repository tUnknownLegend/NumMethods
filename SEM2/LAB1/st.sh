mkdir build
cd build
scan-build cmake -DCMAKE_BUILD_TYPE=Debug ..
scan-build cmake ..
scan-build make
#clear
mkdir ./CMakeFiles/exec
mv LAB4 ./CMakeFiles/exec
cd ./CMakeFiles/exec
./LAB4
cd ../../../
rm -r ./build