#!/bin/bash

cuda_path_inc=/usr/local/cuda/include
cuda_path_lib=/usr/local/cuda/lib64

nvcc -c source/routines.cu -Iinclude -I$cuda_path_inc -DPPT_ENABLE_CUDA_BACKEND
nvcc -o main.exe main.cpp routines.o -Iinclude -I$cuda_path_inc -L$cuda_path_lib -lcudart -DPPT_ENABLE_CUDA_BACKEND
