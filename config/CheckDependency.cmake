##############################################################################
# @file  checkDependency.cmake
# @brief Macro definition file.
#
# Copyright (c) 2013 University of Pennsylvania. All rights reserved.<br />
# See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

macro(check_dependency command NOT_FOUND_TYPE)
	# NOT_FOUND_TYPE paramemter is the keyword for message when the command is not found. It is identical to the keywords for message() function. See: http://www.cmake.org/cmake/help/v2.8.4/cmake.html#command:message
	string(TOUPPER ${command} _COMMAND)
	find_program(${_COMMAND}_PATH ${command})
	if (${_COMMAND}_PATH STREQUAL ${_COMMAND}_PATH-NOTFOUND) 
		if (${NOT_FOUND_TYPE} STREQUAL SEND_ERROR OR ${NOT_FOUND_TYPE} STREQUAL FATAL_ERROR)
			message(${NOT_FOUND_TYPE} "${command} is REQUIRED to run the software however it cannot be found in your enviornment variable PATH. Please install ${command}, or check your PATH variable if you have already installed ${command}.")
		elseif(${NOT_FOUND_TYPE} STREQUAL WARNING)
			message(${NOT_FOUND_TYPE} "${command} is an optional dependency for ${CMAKE_PROJECT_NAME}. You may still use the software properly with appropriate options. For more information, check the help documentation.")
		endif()
	else()
		set(${_COMMAND}_FOUND 1)
	endif ()
endmacro()

