cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(chrp CXX)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)

find_package(Chemharp)

if(NOT ${CHEMHARP_FOUND})
    message(FATAL_ERROR "Could not find Chemharp. Please set CHEMHARP_INCLUDEDIR \
                         and CHEMHARP_LIBDIR to the location of the library.")
endif()

include_directories(${CHEMHARP_INCLUDE_DIRS})

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'release' as none was specified.")
  set(CMAKE_BUILD_TYPE "release" CACHE STRING "Choose the type of build." FORCE)
endif()


include_directories(SYSTEM external/)
add_library(docopt STATIC external/docopt/docopt.cpp)

find_library(CHEMHARP_LIBRARY NAMES chemharp)
find_library(CHEMHARP_LIBRARY NAMES chemharp)

file(GLOB_RECURSE sources src/lib/**.cpp)
add_library(frontend STATIC ${sources})

add_executable(chrp src/main.cpp)
target_link_libraries(chrp frontend docopt ${CHEMHARP_LIBRARIES})


install(TARGETS chrp DESTINATION bin)

# uninstall target
configure_file(
    "${CMAKE_SOURCE_DIR}/cmake/uninstall.in.cmake"
    "${CMAKE_BINARY_DIR}/uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake)
