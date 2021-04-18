#############################################################
#############################################################
# Print all arguments
macro (print_all_args)
  message (STATUS "------------------------------------")
  
  set(list_var "${ARGN}")
  foreach(_currentItem IN LISTS list_var)
    message (STATUS "Adding library source file:   ${_currentItem}")
  endforeach (_currentItem)
  
  message (STATUS "------------------------------------")
endmacro (print_all_args)
#############################################################
#############################################################



###############################################################################
###############################################################################
message(STATUS "This is CMake ${CMAKE_VERSION}")
message(STATUS "")

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

FIND_PACKAGE(deal.II 9.3.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#
# Are all dependencies fulfilled?
#
IF(NOT ((DEAL_II_WITH_PETSC AND NOT DEAL_II_PETSC_WITH_COMPLEX) OR DEAL_II_WITH_TRILINOS) OR NOT DEAL_II_WITH_P4EST) # keep in one line
  MESSAGE(FATAL_ERROR "
Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_PETSC = ON
    DEAL_II_PETSC_WITH_COMPLEX = OFF
    DEAL_II_WITH_P4EST = ON
or
    DEAL_II_WITH_TRILINOS = ON
    DEAL_II_WITH_P4EST = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
    DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
    DEAL_II_PETSC_WITH_COMPLEX = ${DEAL_II_PETSC_WITH_COMPLEX}
    DEAL_II_WITH_P4EST = ${DEAL_II_WITH_P4EST}
    DEAL_II_WITH_TRILINOS = ${DEAL_II_WITH_TRILINOS}
which conflict with the requirements.
One or both of the aforementioned combinations of prerequisites are not met by your installation, but at least one is required for this tutorial step."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug."
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Release."
  )
  
ADD_CUSTOM_TARGET(indent
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./contrib/utilities/indent
  COMMENT "Indenting all C++ source and header files."
  )
  
ADD_CUSTOM_TARGET(documentation
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND ./contrib/utilities/build_docs
  COMMENT "Building doxygen documentation."
  ) 
###############################################################################
###############################################################################


#####################################
#####################################
#
# Now our stuff
#
set(PROJECT_NAME P4EST_INTERFACE)

# The version number
set(P4EST_INTERFACE_VER_MAJOR 0)
set(P4EST_INTERFACE_VER_MINOR 1)

project(${PROJECT_NAME})

# Check for the existence of various optional folders:
if (EXISTS ${CMAKE_SOURCE_DIR}/doc/CMakeLists.txt)
  set(P4EST_INTERFACE_HAVE_DOC_DIRECTORY TRUE)
endif ()

if (EXISTS ${CMAKE_SOURCE_DIR}/test/CMakeLists.txt)
  set (P4EST_INTERFACE_HAVE_TEST_DIRECTORY TRUE)
endif ()


# Change default CMAKE_INSTAL_PREFIX to ${CMAKE_BINARY_DIR}/lib
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "default install path" FORCE )
endif ()

set (P4EST_INTERFACE_LIBRARY P4EST_INTERFACE)


set(P4EST_INTERFACE_INCLUDE_DIR
		${CMAKE_SOURCE_DIR}/include)

set(P4EST_INTERFACE_SRC_DIR 
		${CMAKE_SOURCE_DIR}/source)
add_subdirectory (${P4EST_INTERFACE_SRC_DIR})


# Only add subdirectories for doc if exists
if (P4EST_INTERFACE_HAVE_DOC_DIRECTORY)
  add_subdirectory (${CMAKE_SOURCE_DIR}/doc)
endif ()

# Only add subdirectories for doc if exists
if (P4EST_INTERFACE_HAVE_TEST_DIRECTORY)
  add_subdirectory (${CMAKE_SOURCE_DIR}/test)
  ENABLE_TESTING()
endif ()



#####################################
#####################################