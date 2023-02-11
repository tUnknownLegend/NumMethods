mkdir build
cd build
scan-build cmake -DCMAKE_BUILD_TYPE=Debug ..
scan-build cmake ..
scan-build make
#clear
mkdir ./CMakeFiles/exec
mv LAB2 ./CMakeFiles/exec
cd ./CMakeFiles/exec
./LAB2
cd ../../../
rm -r ./build