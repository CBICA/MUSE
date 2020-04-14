##############################################################################
# @file  Settings.cmake
# @brief Non-default project settings.
#
# This file is included by basis_project_impl(), after it looked for the
# required and optional dependencies and the CMake variables related to the
# project directory structure were defined (see Directories.cmake file in
# @c BINARY_CONFIG_DIR). It is further included before the BasisSettings.cmake
# file.
#
# In particular build options should be added in this file using CMake's
# <a href="http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:option">
# option()</a> command. Further, any common settings related to using a found
# dependency can be set here if the basis_use_package() command was enable
# to import the required configuration of a particular external package.
#
# Copyright (c) 2013 University of Pennsylvania. All rights reserved.<br />
# See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
#
# @ingroup BasisSettings
##############################################################################

set (DART_TESTING_TIMEOUT 3600 CACHE STRING "Maximum time allowed before CTest will kill the test." FORCE)
set (BUILD_TESTING ON CACHE BOOL "Build the testing tree." FORCE)
set (SCHEDULER "NONE" CACHE STRING "The manufacturer of the cluster job scheduler. Enter one of these {NONE, SGE, PBS, MISC}. Default: NONE")
set (CMAKE_C_FLAGS "-w" CACHE STRING "Flags used by the compiler during all build types." FORCE)

if(SCHEDULER)
	string(TOUPPER ${SCHEDULER} SCHEDULER_U)
	if (SCHEDULER_U MATCHES "SGE")
		set (SCHEDULER_SETTING_FILE "SettingsSGE.sh")
	elseif (SCHEDULER_U MATCHES "PBS")
		set (SCHEDULER_SETTING_FILE "SettingsPBS.sh")
	elseif (SCHEDULER_U MATCHES "MISC")
		set (SCHEDULER_SETTING_FILE "SettingsMISC.sh")
	endif()
endif()
