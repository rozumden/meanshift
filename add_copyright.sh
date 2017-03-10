#!/bin/bash
for i in `find . -name "*.cpp" \
              -o -name "*.c" \
              -o -name "*.cc" \
              -o -name "*.h" \
              -o -name "*.hpp"  -type f`; do
	echo $i
    less copyright_cpp.txt | cat - $i > temp && mv temp $i
done
for i in `find . -name "*.m" -type f`; do
	echo $i
    less copyright_matlab.txt | cat - $i > temp && mv temp $i
done