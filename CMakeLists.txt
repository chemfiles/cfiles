cmake_minimum_required(VERSION 2.8.7)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/chemfiles/cmake")

project(cfiles CXX)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)

if(EXISTS "${PROJECT_SOURCE_DIR}/chemfiles/CMakeLists.txt")
    add_subdirectory(chemfiles)
else()
    message(FATAL_ERROR "Missing the chemfiles library. Please run `git submodule update --init`.")
endif()

include(CompilerFlags)

include_directories(chemfiles/include/)

include_directories(SYSTEM external/)
add_library(docopt STATIC external/docopt/docopt.cpp)

file(GLOB_RECURSE sources src/lib/**.cpp)
add_library(frontend STATIC ${sources})

add_executable(cfiles src/main.cpp)
target_link_libraries(cfiles frontend docopt)

target_link_libraries(cfiles chemfiles)

install(TARGETS cfiles DESTINATION bin)
