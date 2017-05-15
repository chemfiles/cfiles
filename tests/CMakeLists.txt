add_library(catch OBJECT external/catch.cpp)

function(unit_test _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_executable(${_name_} ${_file_} $<TARGET_OBJECTS:catch>)
    target_link_libraries(${_name_} libcfiles)
    target_include_directories(${_name_} PRIVATE external)
    add_test(${_name_} ${_name_})
    set_tests_properties(${_name_} PROPERTIES LABELS "cfiles")
endfunction()

file(GLOB all_unit_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

foreach(test_file IN LISTS all_unit_tests)
    unit_test(${test_file})
endforeach()

###############################################################################

find_package(PythonInterp)
if(NOT ${PYTHONINTERP_FOUND})
    message(ERROR "Python is needed to run the cfiles tests")
endif()

function(cfiles_test _file_)
    get_filename_component(_name_ ${_file_} NAME_WE)
    add_test(
        NAME ${_name_}
        COMMAND ${PYTHON_EXECUTABLE} ${_file_}
        WORKING_DIRECTORY ${cfiles_BINARY_DIR}
    )
    set_tests_properties(${_name_} PROPERTIES LABELS "cfiles")
endfunction()

file(GLOB all_python_tests
    ${CMAKE_CURRENT_SOURCE_DIR}/*.py
)

foreach(test_file IN LISTS all_python_tests)
    cfiles_test(${test_file})
endforeach()
