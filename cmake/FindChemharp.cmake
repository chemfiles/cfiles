# - Try to find Chemharp
#
# This module set the following variables:
#  CHEMHARP_FOUND          - Chemharp was found
#  CHEMHARP_INCLUDE_DIRS   - The Chemharp include directories
#  CHEMHARP_LIBRARIES      - The library needed to use Chemharp
#
# You can use the following as hints for the search:
#  CHEMHARP_INCLUDEDIR     - Where to look for the headers
#  CHEMHARP_LIBDIR         - Where to look for the library
#

find_path(
    CHEMHARP_INCLUDE_DIR chemharp.hpp
    HINTS ${CHEMHARP_INCLUDEDIR}
)

find_library(
    CHEMHARP_LIBRARY
    NAMES chemharp
    HINTS ${CHEMHARP_LIBDIR}
)

set(CHEMHARP_LIBRARIES ${CHEMHARP_LIBRARY})
set(CHEMHARP_INCLUDE_DIRS ${CHEMHARP_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(chemharp DEFAULT_MSG CHEMHARP_LIBRARY CHEMHARP_INCLUDE_DIR)

mark_as_advanced(CHEMHARP_INCLUDE_DIR CHEMHARP_LIBRARY)
