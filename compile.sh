#!/bin/bash

g++ -c source/routines.cpp -Iinclude -fopenmp -DPPT_ENABLE_OPENMP_BACKEND
g++ -o main.exe main.cpp routines.o -Iinclude -fopenmp -DPPT_ENABLE_OPENMP_BACKEND

