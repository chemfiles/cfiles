cmake_minimum_required(VERSION 2.8.11)

project(cfiles C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CXX_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones" OFF)

if(EXISTS "${PROJECT_SOURCE_DIR}/chemfiles/CMakeLists.txt")
    add_subdirectory(chemfiles)
else()
    message(FATAL_ERROR "Missing the chemfiles library. Please run `git submodule update --init`.")
endif()

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Git QUIET)
set(CFILES_VERSION 0.1.0)

file(WRITE ${CMAKE_BINARY_DIR}/version.hpp.in
"\#define CFILES_VERSION \"\${CFILES_VERSION}\"\n"
)
file(WRITE ${CMAKE_BINARY_DIR}/version.cmake
"
cmake_policy(SET CMP0012 NEW)
set(CFILES_VERSION ${CFILES_VERSION})
if(${GIT_FOUND})
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
        OUTPUT_VARIABLE GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (NOT \"\${GIT_VERSION}\" STREQUAL \"${CFILES_VERSION}\")
        set(CFILES_VERSION ${CFILES_VERSION}-\${GIT_VERSION})
    endif()
endif()
configure_file(\${SRC} \${DST})
")

add_custom_target(version
    ${CMAKE_COMMAND} -D SRC=${CMAKE_BINARY_DIR}/version.hpp.in
                     -D DST=${CMAKE_BINARY_DIR}/version.hpp
                     -P ${CMAKE_BINARY_DIR}/version.cmake
)


add_library(docopt OBJECT external/docopt/docopt.cpp)
target_include_directories(docopt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external)

add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external)
target_compile_definitions(eigen INTERFACE -DEIGEN_MPL2_ONLY)

add_library(kissfft OBJECT external/kissfft/kiss_fft.c external/kissfft/tools/kiss_fftr.c)
target_include_directories(kissfft PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/kissfft
    ${CMAKE_CURRENT_SOURCE_DIR}/external/kissfft/tools
)

file(GLOB_RECURSE sources ${CMAKE_CURRENT_SOURCE_DIR}/src/**.cpp)
list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

add_library(libcfiles
    ${sources}
    $<TARGET_OBJECTS:docopt>
    $<TARGET_OBJECTS:kissfft>
    $<TARGET_OBJECTS:fmt::object>
)
set_target_properties(libcfiles PROPERTIES OUTPUT_NAME cfiles)

target_include_directories(libcfiles
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_BINARY_DIR}/
    PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/chemfiles/external/fmt/include
        ${CMAKE_CURRENT_SOURCE_DIR}/external/kissfft
        ${CMAKE_CURRENT_SOURCE_DIR}/external/kissfft/tools
)

target_link_libraries(libcfiles eigen chemfiles)

include(CompilerFlags)
try_compile(STD_REGEX_WORKS
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/cmake/check_std_regex.cpp
)

if(NOT ${STD_REGEX_WORKS})
    message(STATUS "Using Boost.Regex instead of std::regex")
    find_package(Boost 1.53 REQUIRED COMPONENTS regex)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(docopt ${Boost_LIBRARIES})
    target_compile_definitions(docopt PRIVATE "-DDOCTOPT_USE_BOOST_REGEX")
endif()

add_dependencies(libcfiles version)

add_executable(cfiles src/main.cpp)
target_link_libraries(cfiles libcfiles)

install(TARGETS cfiles DESTINATION bin)

if(CHFL_CODE_COVERAGE)
    # Code coverage should use gcc
    if(NOT CMAKE_COMPILER_IS_GNUCXX)
        message(FATAL_ERROR "Code coverage can only be used with GCC")
    endif()
    set(CMAKE_BUILD_TYPE "debug" CACHE STRING "Type of build: release|debug|relwithdebuginfo" FORCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
endif()

enable_testing()
add_subdirectory(tests)
