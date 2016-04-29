cmake_minimum_required(VERSION 2.8.11)
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

file(GLOB sources src/*.cpp)
add_executable(cfiles ${sources})
target_link_libraries(cfiles docopt chemfiles)

if("${CMAKE_COMPILER_IS_GNUCXX}" AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.9)
    message(STATUS "Using Boost.Regex instead of std::regex")
    find_package(Boost REQUIRED COMPONENTS Regex)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(cfiles ${Boost_LIBRARIES})
    add_definitions("-DDOCOPT_USES_BOOST_REGEX")
endif()

install(TARGETS cfiles DESTINATION bin)

enable_testing()
add_subdirectory(tests)

if(CODE_COVERAGE)
    # Code coverage should use gcc
    if(NOT CMAKE_COMPILER_IS_GNUCXX)
        message(FATAL_ERROR "Code coverage can only be used with GCC")
    endif()
    set(CMAKE_BUILD_TYPE "debug" CACHE STRING "Type of build: release|debug|relwithdebuginfo" FORCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
endif()
