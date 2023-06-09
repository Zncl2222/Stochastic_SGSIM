#!/bin/bash

is_execute=true

# Parse command-line arguments
while getopts ":s" opt; do
    case $opt in
        s)
            is_execute=false
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            ;;
    esac
done

mkdir cmake_build
cd cmake_build

if [ "$is_execute" = true ]; then
    echo "Compile the code into an executable file !"
    cmake -DIS_EXECUTE=ON ..
else
    echo "Compile the code into a shared library (dynamic library) !"
    cmake -DIS_EXECUTE=OFF ..
fi

make
make test
make memcheck
