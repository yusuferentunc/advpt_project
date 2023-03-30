#!/bin/bash -l
#module load cmake
cmake .

echo "The project should be built now..."

make clean
make