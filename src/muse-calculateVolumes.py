#!/usr/bin/env python

import os, sys

EXEC_NAME = "muse-calculateVolumes"

############################ FUNCTIONS ############################

### Version info
def print_version():
	ver = "$Revision: 391 $"
	print("@PROJECT_VERSION@ ( rev" + ver[:-1] + ")")
	sys.exit(0)

# Usage info
def usage():
    """usage information"""
    print(r"""
  %(EXEC)s--
    Calculate volumes of a given image

  Usage: %(EXEC)s [OPTIONS]

  OPTIONS:
  
  Basic:
    [-i --input     ]	< file	>	Input image
    [-M --mask	    ]	< file >	Mask image in the --input image space (optional)
    [-s --subid     ]	< ID	>	Subject ID (optional)
    [-o --output    ]	< file	>	Output csv file (optional)
    [-n --nonzero   ]			Calculate volumes for the non-zero intensities (optional)
    [-I --intensity ]	< int	>	Calculate volumes for the given intensity (optional)

  Derived:
    [-d --derived   ]			Calculate derived volumes as well (optional)
    [-v --ICV	    ]	< file >	File to calculate the ICV from (optional, valid only if calculating derived volumes)
    [-m --map	    ]	< file  >	Mapping for all derived ROIs (optional, valid only if calculating derived volumes)
    					default: '@DATA_DIR@/List/MUSE_DerivedROIs_Mappings.csv'
    [-V --version   ]			Version information

  """ % {'EXEC':EXEC_NAME})

def signal_handler(signal, frame):
        print('You pressed Ctrl+C!')
        pool.terminate()
        sys.exit(0)

def calculateVolume(roi):
	if not maskfile:
		return voxvol * np.count_nonzero( img[ img==roi ] )
	else:
		return voxvol * np.count_nonzero( img[ mask==roi ] )

############################ MAIN ############################

# Check the number of arguments
if len(sys.argv) < 2:
	usage()
	sys.exit(0)

# Checking if the first argument is -V
if sys.argv[1] in ["-V", "--version"]:
	print_version()

### Import some more modules
import nibabel as nib
import numpy as np
import signal, getopt

### Default parameters
maskfile 	= None
ID 		= None
csvfile 	= None
icvfile 	= None
unique 		= int(1)
nonzero 	= None
intensity 	= None
derived 	= int(0)
ROIs 		= []
ROIlist 	= []
Vollist 	= []
mapcsv 		= '@DATA_DIR@/List/MUSE_DerivedROIs_Mappings.csv'

### Read command line args
try:
	opts, files = getopt.getopt(sys.argv[1:], "i:s:o:nI:dm:v:M:V",
	  ["input", "subid", "output", "nonzero", "intensity", "derived", "map", "ICV", "mask", "version"])

except getopt.GetoptError as err:
	usage()
	print("ERROR!", err)

for o, a in opts:
	if o in ["-i", "--input"]:
		datafile = a
	elif o in ["-M", "--mask"]:
		maskfile = a
	elif o in ["-s", "--subid"]:
		ID = a
	elif o in ["-o", "--output"]:
		csvfile = a
	elif o in ["-n", "--nonzero"]:
		nonzero = int(1)
	elif o in ["-I", "--intensity"]:
		intensity = int(a)
	elif o in ["-d", "--derived"]:
		derived = int(1)
	elif o in ["-m", "--map"]:
		mapcsv = a
	elif o in ["-v", "--ICV"]:
		icvfile = a
	elif o in ["-V", "--version"]:
		print_version
	else:
		usage()
		sys.exit(0)

### Sanity check on the arguments
if not datafile:
	print("ERROR: Input file not provided!!!")
	sys.exit(0) 

if derived == 1 and not mapcsv:
	print("ERROR: Map csv file not provided!!!")
	sys.exit(0) 

#if derived == 1 and not icvfile:
#	print "ERROR: File to calculate ICV not provided!!!"
#	sys.exit(0) 

### Read the input image
data = nib.load(os.path.join(datafile))
img = data.get_data()
hdr = data.get_header()

### Read the mask file
if maskfile:
	mask = nib.load( os.path.join(maskfile) ).get_data()

### Retaining only non-zero values from image
if not maskfile:
	img = img[ img>0 ]
else:
	maskInd = np.where( mask > 0 )
	mask = mask[ maskInd ]
	img = img[ maskInd ]

### Get voxel dimensions
voxdims = hdr.structarr["pixdim"]
voxvol = float(voxdims[1]*voxdims[2]*voxdims[3])

### Reading the data from icv file and getting voxel dims
if derived == 1 and icvfile:
	icvdata = nib.load(os.path.join(icvfile))
	icvimg = icvdata.get_data()
	icvhdr = icvdata.get_header()

	icvvoxdims = icvhdr.structarr["pixdim"]
	icvvoxvol = float(icvvoxdims[1]*icvvoxdims[2]*icvvoxdims[3])
	
	
### Decide on the ROIs to calculate volume for
if intensity:
	ROIs = [ intensity ]
elif nonzero and nonzero == 1:
	img[ img>0 ] = 1
	ROIs = [ int(1) ]
elif unique and unique == 1:
	if not maskfile:
		ROIs = np.unique(img[np.nonzero(img)])
	else:
		ROIs = np.unique(mask[np.nonzero(mask)])

### Calculate volumes
Volumes = []
Volumes = np.array( list(map( calculateVolume, (ROIs) )) )

DerivedROIs = []
DerivedVols = []

### Calculate derived volumes if requested
if derived == 1:
	import csv
	
#	with open(mapcsv, 'rb') as mapcsvfile:
	with open(mapcsv) as mapcsvfile:
		reader = csv.reader(mapcsvfile, delimiter=',')
		
		# Read each line in the csv map files
		for row in reader:			
			# Calculate ICV
			if row[0] == '702':
				if icvfile:
					# Append the ROI number to the list
					DerivedROIs.append(row[0])
					vol = icvvoxvol*np.count_nonzero(icvimg[ icvimg>0 ])
					DerivedVols.append(vol)

			# Calculate volumes for other derived ROIs
			else:
				# Append the ROI number to the list
				DerivedROIs.append(row[0])

				vol = 0
				for roi in range( 2,len(row) ):
					# Check if the roi exists in the list. If it does, calculate, else just say 0
					if ROIs[ ROIs == int(row[ roi ]) ]:
						vol += Volumes[ ROIs == int(row[ roi ]) ][0]
					
				DerivedVols.append(vol)

### Decide what to output
if ID:
	ROIlist = ['ID']
	Vollist = [ID]

if derived == 1:
	ROIlist.extend(DerivedROIs)
	Vollist.extend(DerivedVols)
else:
	ROIlist.extend(ROIs)
	Vollist.extend(Volumes)

### Decide where to output
if csvfile:
	with open(csvfile, 'w+') as f:
		f.write(str(','.join(map(str, ROIlist))))
		f.write(str('\n'))
		f.write(str(','.join(map(str, Vollist))))
		f.write(str('\n'))
else:
	print(','.join(map(str, ROIlist)))
	print(','.join(map(str, Vollist)))
