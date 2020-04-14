##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of software.
#
# Copyright (c) 2013 University of Pennsylvania. All rights reserved.<br />
# See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

# ============================================================================
# library target(s)
# ============================================================================

# Add library target for each library using basis_add_library().
#
# This command can not only be used to build libraries from C/C++, but also source
# code written in other programming languages such as Java, Python, Perl,
# MATLAB, and Bash. Note that here we consider modules written in a scripting
# language which are no executables but to be included by other scripts written
# in the particular language as libraries.
#
# Note: Public modules written in a scripting language such as Python, Perl,
#       MATLAB, or Bash which are intended for use by other packages should
#       be placed in the lib/[<lang>/]MUSE/ directory,
#       where <lang> is the language name in lowercase and is optional.
#       BASIS will automatically add a library build target for these modules.



# ============================================================================
# executable target(s)
# ============================================================================

# Add executable target for each executable program using basis_add_executable().
#
# This command can not only be used to build executables from C/C++, but also
# source code written in other programming languages such as Java, Python, Perl,
# MATLAB, and Bash.

#basis_add_executable (muse.sh NO_BASIS_UTILITIES)
configure_file( muse.sh.in ${CMAKE_BINARY_DIR}/src/muse.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION bin 
	RENAME muse )

#basis_add_executable (muse-calculateVolumes.py NO_BASIS_UTILITIES)
configure_file( muse-calculateVolumes.py ${CMAKE_BINARY_DIR}/src/muse-calculateVolumes.py @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-calculateVolumes.py 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION bin 
	RENAME muse-calculateVolumes )


#basis_add_executable (muse-functions.sh LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-functions.sh.in ${CMAKE_BINARY_DIR}/src/muse-functions.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-functions.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-functions )
#basis_add_executable (muse-dramms-registrations.sh LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-dramms-registrations.sh.in ${CMAKE_BINARY_DIR}/src/muse-dramms-registrations.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-dramms-registrations.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-dramms-registrations )
#basis_add_executable (muse-ants-registrations.sh LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-ants-registrations.sh.in ${CMAKE_BINARY_DIR}/src/muse-ants-registrations.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-ants-registrations.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-ants-registrations )
#basis_add_executable (muse-calculateSimilarityMap.sh LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-calculateSimilarityMap.sh.in ${CMAKE_BINARY_DIR}/src/muse-calculateSimilarityMap.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-calculateSimilarityMap.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-calculateSimilarityMap )
#basis_add_executable (muse-labelFusion.sh LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-labelFusion.sh.in ${CMAKE_BINARY_DIR}/src/muse-labelFusion.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-labelFusion.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-labelFusion )
#basis_add_executable (muse-combine1ROI.sh LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-combine1ROI.sh.in ${CMAKE_BINARY_DIR}/src/muse-combine1ROI.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-combine1ROI.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-combine1ROI )
#basis_add_executable (muse-combineAll.sh LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-combineAll.sh.in ${CMAKE_BINARY_DIR}/src/muse-combineAll.sh.in @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-combineAll.sh.in 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-combineAll )
#basis_add_executable (muse-combineRoiMapsIter.py LIBEXEC NO_BASIS_UTILITIES)
configure_file( muse-combineRoiMapsIter.py ${CMAKE_BINARY_DIR}/src/muse-combineRoiMapsIter.py @ONLY )
install( FILES ${CMAKE_BINARY_DIR}/src/muse-combineRoiMapsIter.py 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib 
	RENAME muse-combineRoiMapsIter )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm")

#basis_add_library(simlib STATIC
#	IOZhan.c
#	cres.c
#	imageZHAN.c
#	matrixZHAN.c
#	mvcd.c
#	watershed.c)
add_library(simlib STATIC
	IOZhan.c
	cres.c
	imageZHAN.c
	matrixZHAN.c
	mvcd.c
	watershed.c )


install( TARGETS simlib
    	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
    	DESTINATION lib )

#basis_add_executable (calculateSimilarityMap.c LIBEXEC)
add_executable( calculateSimilarityMap calculateSimilarityMap.c )

#basis_target_link_libraries(calculateSimilarityMap simlib)
target_link_libraries( calculateSimilarityMap simlib )

install( FILES ${CMAKE_BINARY_DIR}/src/calculateSimilarityMap 
	PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ 
	DESTINATION lib )

if(SCHEDULER)
	add_subdirectory(schedulerSettings)
endif()

