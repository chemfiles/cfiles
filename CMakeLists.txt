cmake_minimum_required(VERSION 2.8.12)

project(chrp CXX)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE sources src/lib/**.cpp)
add_library(frontend STATIC ${sources})

add_executable(chrp src/main.cpp)
target_link_libraries(chrp frontend chemharp)
