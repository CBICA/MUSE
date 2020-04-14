##############################################################################
# @file  CMakeLists.txt
# @brief Root build configuration file.
#
# Copyright (c) 2013 University of Pennsylvania. All rights reserved.<br />
# See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

# ----------------------------------------------------------------------------
# minimum required CMake version
cmake_minimum_required( VERSION 2.8.4 )

project( MUSE )


# ============================================================================
# FIND PACKAGES
# ============================================================================

# include MUSE macros
include(${CMAKE_SOURCE_DIR}/config/CheckDependency.cmake)

list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/config" )

find_package( DRAMMS REQUIRED )
find_package( ANTS REQUIRED )
find_package( AFNITools REQUIRED )
find_package( MICO REQUIRED )
find_package( FSLTools REQUIRED )


# ============================================================================
# SET VARIABLES FOR CMAKE
# ============================================================================

#set (CMAKE_INSTALL_PREFIX "" CACHE INTERNAL "Installation prefix." FORCE)
set (SCHEDULER "NONE" CACHE STRING "The manufacturer of the cluster job scheduler. Enter one of these {NONE, SGE, PBS, MISC}. Default: NONE")
set (CMAKE_C_FLAGS "-w" CACHE STRING "Flags used by the compiler during all build types." FORCE)

set (INSTALL_LIBRARY_DIR lib)
set (INSTALL_LIBEXEC_DIR lib)
set (INSTALL_DATA_DIR data)
set (INSTALL_BASH_LIBRARY_DIR lib/bash)
set (PROJECT_NAMESPACE_BASH muse)
set (INSTALL_DOC_DIR doc)

set (DATA_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_DATA_DIR}")
set (LIBRARY_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBRARY_DIR}")
set (LIBEXEC_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBEXEC_DIR}")
set (BASH_LIBRARY_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_BASH_LIBRARY_DIR}/${PROJECT_NAMESPACE_BASH}")

if(SCHEDULER)
	string(TOUPPER ${SCHEDULER} SCHEDULER_U)
	if (SCHEDULER_U MATCHES "SGE")
		find_package( QSUB REQUIRED )
		set (SCHEDULER_SETTING_FILE "SettingsSGE.sh")
	elseif (SCHEDULER_U MATCHES "PBS")
		set (SCHEDULER_SETTING_FILE "SettingsPBS.sh")
	elseif (SCHEDULER_U MATCHES "MISC")
		set (SCHEDULER_SETTING_FILE "SettingsMISC.sh")
	endif()
endif()

# ============================================================================
# ADD SUBDIRECTORIES
# ============================================================================

add_subdirectory(src)
add_subdirectory(data)

#add_subdirectory(test)


# ============================================================================
# copy the manual pdf from source to installdir during install
# ============================================================================

set (MANUAL_OUTPUT_NAME "${PROJECT_NAME}_Software_Manual")
set (PROJECT_DOC_DIR "${PROJECT_SOURCE_DIR}/doc")

if (EXISTS "${PROJECT_DOC_DIR}/${MANUAL_OUTPUT_NAME}.pdf")
  install (
    FILES            ${PROJECT_DOC_DIR}/${MANUAL_OUTPUT_NAME}.pdf
    DESTINATION      ${INSTALL_DOC_DIR}/
  )

  install (
    FILES            ${PROJECT_SOURCE_DIR}/COPYING.txt
    DESTINATION      ${INSTALL_DOC_DIR}/
  )
endif ()


## ============================================================================
## print all variables
## ============================================================================

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
