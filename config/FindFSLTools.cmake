##############################################################################
# @file  FindFSLTools.cmake
# @brief Try to find FSL tools 
#
# Copyright (c) 2013 University of Pennsylvania. All rights reserved.<br />
# See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

check_dependency(fslmerge SEND_ERROR)
check_dependency(fslsplit SEND_ERROR)

set(ALLFOUND FSLMERGE_FOUND AND FSLSPLIT_FOUND)
if (${ALLFOUND})
	set(FSLTools_FOUND 1)
endif(${ALLFOUND})