#
#  Try to find Mini-XML library
#
#  The module defines:
#       MXML_FOUND        - system has Mini-XML library
#       MXML_LIBRARIES    - libraries (mxml or mxml1)
#       MXML_INCLUDE_DIR  - path to header file
#

include(FindPackageHandleStandardArgs)

find_package(Threads REQUIRED)

if (UNIX) # assume on the UNIX platform the Mini-XML library was installed in standard way
    find_path (MXML_INCLUDE_DIR mxml.h)
    find_library (MXML_LIBRARY NAMES mxml PATH_SUFFIXES lib lib64)
elseif (WIN32) # user may set MXML_INSTALL_DIR variable to detect header and library files
    set(MXML_INSTALL_DIR "" CACHE STRING "MXML install dir")
    set(MXML_INSTALL_DIR_INTERNAL "" CACHE STRING "MXML install dir")
    if(NOT "${MXML_INSTALL_DIR}" STREQUAL "${MXML_INSTALL_DIR_INTERNAL}") # MXML_INSTALL_DIR is given in command-line
        unset(MXML_LIBRARY CACHE)
        unset(MXML_INCLUDE_DIR CACHE)
    endif()
    find_path (MXML_INCLUDE_DIR NAMES mxml.h PATHS ${MXML_INSTALL_DIR})
    find_library (MXML_LIBRARY NAMES mxml mxml1 PATHS ${MXML_INSTALL_DIR})
endif ()

    find_package_handle_standard_args(MXML  DEFAULT_MSG
                                      MXML_LIBRARY MXML_INCLUDE_DIR)
    mark_as_advanced(MXML_INCLUDE_DIR MXML_LIBRARY)
    set(MXML_LIBRARIES ${MXML_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
