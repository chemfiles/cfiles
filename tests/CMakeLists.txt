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

file(GLOB all_test_files
    ${CMAKE_CURRENT_SOURCE_DIR}/*.py
)

foreach(test_file IN LISTS all_test_files)
    cfiles_test(${test_file})
endforeach(test_file)
