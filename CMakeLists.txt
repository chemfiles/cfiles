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
include_directories(src/)
include_directories(SYSTEM external/)

file(GLOB_RECURSE sources src/**.cpp)
add_executable(cfiles ${sources})

try_compile(STD_REGEX_WORKS
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/check_std_regex.cpp
)

if(NOT ${STD_REGEX_WORKS})
    message(STATUS "Using Boost.Regex instead of std::regex")
    find_package(Boost 1.53 REQUIRED COMPONENTS regex)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(cfiles ${Boost_LIBRARIES})
    add_definitions("-DDOCTOPT_USE_BOOST_REGEX")
endif()

add_library(docopt STATIC external/docopt/docopt.cpp)
target_link_libraries(cfiles docopt chemfiles)

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
