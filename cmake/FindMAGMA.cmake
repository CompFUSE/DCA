# - Find the MAGMA library
#
# Usage:
#   FIND_PACKAGE(MAGMA [REQUIRED] [QUIET] )
#
# It sets the following variables:
#   MAGMA_FOUND               ... true if MAGMA is found on the system
#   MAGMA_LIBRARIES           ... full path to MAGMA library
#   MAGMA_INCLUDE_DIRS        ... MAGMA include directory
#

include(FindPackageMessage)

SET(MAGMA_ROOT CACHE STRING
        "Root directory for MAGMA implementation")

# Check if we can use PkgConfig
FIND_PACKAGE(PkgConfig)

# Determine from PKG
IF(PKG_CONFIG_FOUND AND NOT MAGMA_ROOT)
    PKG_CHECK_MODULES( PC_MAGMA QUIET "magma")
ENDIF()

IF(PC_MAGMA_FOUND)
    FOREACH(PC_LIB ${PC_MAGMA_LIBRARIES})
        FIND_LIBRARY(${PC_LIB}_LIBRARY NAMES ${PC_LIB} HINTS ${PC_MAGMA_LIBRARY_DIRS} )
        IF (NOT ${PC_LIB}_LIBRARY)
            MESSAGE(FATAL_ERROR "Something is wrong in your pkg-config file - lib ${PC_LIB} not found in ${PC_MAGMA_LIBRARY_DIRS}")
        ENDIF (NOT ${PC_LIB}_LIBRARY)
        LIST(APPEND MAGMA_LIB ${${PC_LIB}_LIBRARY})
    ENDFOREACH(PC_LIB)

    FIND_PATH(
            MAGMA_INCLUDE_DIRS
            NAMES "magma.h"
            PATHS
                ${PC_MAGMA_INCLUDE_DIRS}
                ${INCLUDE_INSTALL_DIR}
                /usr/include
                /usr/local/include
                /sw/include
                /opt/local/include
            DOC "MAGMA Include Directory"
    )

    FIND_PACKAGE_HANDLE_STANDARD_ARGS(MAGMA DEFAULT_MSG MAGMA_LIB)
    MARK_AS_ADVANCED(MAGMA_INCLUDE_DIRS MAGMA_LIB)

ELSE(PC_MAGMA_FOUND)
    IF(MAGMA_DIR AND NOT MAGMA_ROOT)
        set(MAGMA_ROOT "${MAGMA_DIR}")
    ENDIF()

    IF(MAGMA_ROOT)
        #find libs
        FIND_LIBRARY(
                MAGMA_LIB
                NAMES "magma" "MAGMA"
                PATHS
                    ${MAGMA_ROOT}
                PATH_SUFFIXES "lib" "lib64" "lib/ia32" "lib/intel64"
                DOC "MAGMA Library"
                NO_DEFAULT_PATH
        )

        FIND_PATH(
                MAGMA_INCLUDE_DIRS
                NAMES "magma.h"
                PATHS
                    ${MAGMA_ROOT}
                HINTS
                    ENV OLCF_MAGMA_ROOT
                    ENV MAGMADIR
                    ENV MAGMA_ROOT
                    ENV MAGMA_ROOT_DIR
                PATH_SUFFIXES "include"
                DOC "MAGMA Include Directory"
                NO_DEFAULT_PATH
        )
    ELSE()
        FIND_LIBRARY(
                MAGMA_LIB
                NAMES "magma"
                PATHS
                    ${PC_MAGMA_LIBRARY_DIRS}
                    ${LIB_INSTALL_DIR}
                    /usr/lib64
                    /usr/lib
                    /usr/local/lib64
                    /usr/local/lib
                    /sw/lib
                    /opt/local/lib
                HINTS
                    ENV OLCF_MAGMA_ROOT
                    ENV MAGMADIR
                    ENV MAGMA_DIR
                    ENV MAGMA_ROOT
                    ENV MAGMA_ROOT_DIR
                PATH_SUFFIXES "lib" "lib64" "lib/ia32" "lib/intel64"
                DOC "MAGMA Library"
        )

        FIND_PATH(
                MAGMA_INCLUDE_DIRS
                NAMES "magma.h"
                PATHS
                    ${PC_MAGMA_INCLUDE_DIRS}
                    ${INCLUDE_INSTALL_DIR}
                    /usr/include
                    /usr/local/include
                    /sw/include
                    /opt/local/include
                HINTS
                    ENV OLCF_MAGMA_ROOT
                    ENV MAGMADIR
                    ENV MAGMA_DIR
                    ENV MAGMA_ROOT
                    ENV MAGMA_ROOT_DIR
                PATH_SUFFIXES
                    "include"
                    "lapacke"
                    DOC "MAGMA Include Directory"
        )
    ENDIF(MAGMA_ROOT)
ENDIF(PC_MAGMA_FOUND)

IF(PC_MAGMA_FOUND OR (MAGMA_LIB))
    SET(MAGMA_LIBRARIES ${MAGMA_LIB})
ENDIF()
IF(MAGMA_INCLUDE_DIRS)
    SET(MAGMA_INCLUDE_DIR ${MAGMA_INCLUDE_DIRS})
ENDIF()

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MAGMA DEFAULT_MSG
        MAGMA_INCLUDE_DIRS MAGMA_LIBRARIES)

FIND_PACKAGE_MESSAGE(MAGMA "Found MAGMA: ${MAGMA_LIBRARIES}"
     "[${MAGMA_INCLUDE_DIRS}]")

MARK_AS_ADVANCED(
        MAGMA_ROOT
        MAGMA_INCLUDE_DIRS
        MAGMA_INCLUDE_DIR
        MAGMA_LIBRARIES
        MAGMA_LIB)
