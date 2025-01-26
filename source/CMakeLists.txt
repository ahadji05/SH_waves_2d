
if(NOT BACKEND_SELECTED)
message(FATAL_ERROR "NO BACKEND HAS BEEN SUCCESSFULLY SELECTED!")
endif()

if(USE_OPENMP)
set(source_list routines.cpp)
elseif(USE_CUDA)
set(source_list routines.cu)
else()
message("FAILED TO SET SOURCE FILES!!!")
endif()

target_sources(main.exe PRIVATE ${source_list})
