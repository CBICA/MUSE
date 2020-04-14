#!/usr/bin/env python
#
# @file  muse_combineRoiMapsIter.py
# @brief Combine roi probability maps for a single subject
#
# Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
# See http://www.cbica.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
#



#Usage 
# ############################################ #
# muse_combineRoiMapsIter.py /Path/To/Input/List.txt /Path/To/Destination/outImgName
################################################
# will read roi files listed in 'List.txt' file
#The list file must have full paths to the files

import nibabel as nib
import numpy as np
import sys
import re
import time

print(str(sys.argv))

InputList=str(sys.argv[1])
DestFile=str(sys.argv[2])

### Sanity check on the arguments
if not InputList or not DestFile:
	print("ERROR: Required input options not provided!!!")
	sys.exit(0) 

### Printing input arguments
print('\n\n')
print('Subject Input List :', InputList)
print('Destination File :', DestFile)
print('\n\n')

### Reading input file first line
f=open(InputList)
fline = f.readline()
f.close()

### Extract roi no
match=re.search('([\w.-]+)ROI_(\d+)_([\w.-]+)', fline)
if match:
	rnos = match.group(2)
	rno = int(rnos)
else:
	print('ERROR: No ROI_{roino} in file name !')
	exit(1)

### Read img, vectorize
img = nib.load(str.rstrip(fline))
a=img.get_data()
b=np.reshape(a,-1)
isize = a.shape
vsize = b.shape

### Set index of voxels belonging to that roi, set also max values
imgMAX = b
imgIND = np.zeros(vsize)
imgIND[b>0] = rno

### Reading input file list
f=open(InputList)
lines = f.readlines()
f.close()
ctr=1

### Combine roi images
for line in lines:
        
	print(line)
	
	### Extract roi no
	match=re.search('([\w.-]+)ROI_(\d+)_([\w.-]+)', line)
	if match:
		rnos = match.group(2)
		rno = int(rnos)
	else:
		print('ERROR: No ROI_{roino} in file name !')
		exit(1)
	
	### Read img, vectorize
	img = nib.load(str.rstrip(line))
	a=img.get_data()
	b=np.reshape(a,-1)

	### Set index of voxels belonging to that roi, set also max values
	imgIND.put((b>imgMAX).nonzero(), rno)
	imgMAX = np.maximum(b,imgMAX)


### Write out img
imgINDM = np.reshape(imgIND,isize)

aff = img.get_affine()
hdr = img.get_header()
#hdr.set_data_dtype(np.int16)
img2 = nib.Nifti1Image(imgINDM, aff, hdr)
img2.to_filename(DestFile);
