##############################################################################
# @file  FindAFNITools.cmake
# @brief Try to find AFNI tools 
#
# Copyright (c) 2013 University of Pennsylvania. All rights reserved.<br />
# See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

check_dependency(3dcalc SEND_ERROR)
check_dependency(3dMean SEND_ERROR)
check_dependency(3dTsort SEND_ERROR)
check_dependency(3dBrickStat SEND_ERROR)
check_dependency(3dclust SEND_ERROR)
check_dependency(nifti1_test SEND_ERROR)

set(ALLFOUND 3DCALC_FOUND AND 3DMEAN_FOUND AND 3DTSORT_FOUND AND 3DBRICKSTAT_FOUND AND 3DCLUST_FOUND AND NIFTI1_TEST_FOUND)
if (${ALLFOUND})
	set(AFNITools_FOUND 1)
endif(${ALLFOUND})
