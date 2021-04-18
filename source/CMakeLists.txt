###############################################################################
###############################################################################
#
# Include directory for sources
#
include_directories(${P4EST_INTERFACE_INCLUDE_DIR})
###############################################################################
###############################################################################


###############################################################################
###############################################################################
#
# Name all sources
#
#file(GLOB_RECURSE CASIS_TARGET_LIB_SRC  "coefficients/*.cc") # source files
#set(CASIS_TARGET_SRC ${CASIS_TARGET_SRC})

set(UTILITIES_SRC
	RandomNumber.cpp
	Utilities.cpp)

set(P4EST_SEARCH_SRC
	PartitionSearch.inst.cpp)

#
# Summarize what we picked up
#
print_all_args (
  ${CASIS_TARGET_SRC}
 )
###############################################################################
###############################################################################


###############################################################################
###############################################################################
#
# Compile and link the sources as SHARED
#
add_library (UTILITIES SHARED ${UTILITIES_SRC})
DEAL_II_SETUP_TARGET(UTILITIES)

add_library (P4EST_SEARCH SHARED ${P4EST_SEARCH_SRC})
DEAL_II_SETUP_TARGET(P4EST_SEARCH)
TARGET_LINK_LIBRARIES(P4EST_SEARCH
					UTILITIES)


#
# Compile the executables
#
add_executable(p4est_search_partition "p4est_search_partition.cxx")
DEAL_II_SETUP_TARGET(p4est_search_partition)
TARGET_LINK_LIBRARIES(p4est_search_partition 	
					P4EST_SEARCH
					UTILITIES)
###############################################################################
###############################################################################
