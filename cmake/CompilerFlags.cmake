include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

set(CMAKE_REQUIRED_QUIET YES)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Manually check for some flags, as some versions of CMake do not support
# `CMAKE_CXX_STANDARD`
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXXOX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    if(MSVC)
        if(MSVC_VERSION LESS 1900)
            message(SEND_ERROR "MSVC < 14.0 is not supported. Please update your compiler or use mingw")
        endif()
    else()
        message(SEND_ERROR "The ${CMAKE_CXX_COMPILER} compiler lacks C++11 support. Use another compiler.")
    endif()
endif()

if(COMPILER_HAS_HIDDEN_VISIBILITY)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fvisibility=hidden")
endif()

if(MSVC)
    add_definitions("/D NOMINMAX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
    set(CMAKE_SHARED_LINKER_FLAGS "/SUBSYSTEM:CONSOLE")
endif()
