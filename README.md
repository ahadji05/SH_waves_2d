# SH_waves_2d: Modelling of shear-horizontal waves in 2 dimensions using a staggered-grid finite-difference method

Below you can find the instructions to run the demo and visualize the output snapshot using *Python matplotlib*.

## 1) running this demo with OpenMP and visualizing the output snapshot is as simple as:
  ./compile.sh
 
  export OMP_NUM_THREADS=#number_of_threads_to_use
  
  ./run.sh

## 2) to run this demo with CUDA:
./compile_cuda.sh

./run.sh

**NOTE:** cuda installation is expected to be in path: ***/usr/local/cuda***, if not, edit the script *compile_cuda.sh* according to the path in your environment.

## 3) to clean up:

./clean.sh

## 4) Bibliography
The modelling algorithm is developed following the papers attached in directory: */bibliography* .
