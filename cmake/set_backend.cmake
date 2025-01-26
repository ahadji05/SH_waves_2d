
if(USE_CUDA)
enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(main.exe PUBLIC CUDA::cudart CUDA::cufft)
endif()
