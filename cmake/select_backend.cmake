
if(USE_OPENMP)
set(BACKEND_SELECTED ON)
target_compile_definitions(main.exe PRIVATE PPT_ENABLE_OPENMP_BACKEND)
elseif(USE_CUDA)
set(BACKEND_SELECTED ON)
target_compile_definitions(main.exe PRIVATE PPT_ENABLE_CUDA_BACKEND)
else()
set(BACKEND_SELECTED OFF)
message(FATAL_ERROR "NO BACKEND WAS SPECIFIED! CHOOSE AMONG: \n -DUSE_OPENMP=ON \n -DUSE_CUDA=ON")
endif()

include(cmake/set_backend.cmake)
