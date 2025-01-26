cmake_minimum_required(VERSION 3.22)

project(SH_ELASTIC_WAVES_2D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(main.exe main.cpp)

include(cmake/select_backend.cmake)

target_include_directories(main.exe PRIVATE include)

add_subdirectory(source)
