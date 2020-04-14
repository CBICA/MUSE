#!/bin/bash

setSubmitCMD()
{
	# Getting the arguments 
	_memory=$1
	_threads=$2
	shift 2
	_holdids="$*"
	
	# Getting the options for PBS scheduler into variables
	# This is done so that the user can adapt these options to match them with his own scheduler
	_OptExec="qsub"						# Executable name for job scheduler
	_OptPriority="-p"					# Option for setting priority
	_OptMemory="-l mem=${_memory}gb"			# String for requesting memory
	_OptTerse=""						# Option for requesting output of a job id only
	_OptJoin="-j oe"					# String for joining the output files
	_OptOutput="-o ${log}\$PBS_JOBNAME-\$PBS_JOBID.log"	# String for specifying output filename
	_OptHoldID="-W depend=afterany:"			# String for requesting a job dependency
	_OptThreads="-l ppn=${_threads}"			# String for requesting multiple processors on a node
	
	# Setting up the basic command and then the other options will be added onto it
	submitCMD="${_OptExec} ${_OptTerse} ${_OptJoin} ${_OptOutput} ${_OptMemory}"
	
	
	# check if priority variable
	if [ -n "${priority}" ]
	then 
		submitCMD="${submitCMD} ${_OptPriority} ${priority}";
	fi

	# check if more than one threads are being requested
	if [ ${_threads} -gt 1 ]
	then 
		submitCMD="${submitCMD} ${_OptThreads}";
	fi


	# check if the user provided a hold jobid
	# Modify this part of the code depending on how the job ids need to be appended
	if [ -n "${_holdids}" ]
	then 
		_holdidsString=`echo ${_holdids} | sed 's/ /:/g'`
		submitCMD="${submitCMD} ${_OptHoldID}${_holdidsString}"
	fi

	echo $submitCMD
}

setDeleteJob()
{
	# Getting the job IDs
	_jobid="$*"
	
	# Getting the options for PBS scheduler into variables
	# This is done so that the user can adapt these options to match them with his own scheduler
	_OptJobDeleteExec="qdel"					# Executable name for deleting jobs
	
	echo ${_OptJobDeleteExec}
}
