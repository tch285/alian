#!/bin/bash

./test_analysis.py -e 1000 	-o run2.root run2_tstruct.yaml run2_file_list.txt --run 2
./test_analysis.py -e 10000 -o run3.root run3_tstruct.yaml run3_file_list.txt --run 3

