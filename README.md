# SH_waves_2d: Modelling of shear-horizontal waves in 2 dimensions using a staggered-grid finite-difference method

Below you can find the instructions to run the demo and visualize the output snapshot using *Python matplotlib*.

## 1) running this demo with OpenMP and visualizing the output snapshot is as simple as:
  mkdir build
  cd build
  cmake .. -DUSE_OPENMP=ON
  export OMP_NUM_THREADS=#nthreads_to_use
  ./../run.sh

## 2) to run this demo with CUDA:
  mkdir build
  cd build
  cmake .. -DUSE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=#cuda_arch
  ./../run.sh

## 3) Bibliography
The modelling algorithm is developed following the papers attached in directory: */bibliography* .
