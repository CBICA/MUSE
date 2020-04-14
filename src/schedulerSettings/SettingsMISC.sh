#!/bin/bash

setSubmitCMD()
{
	# Getting the arguments 
	_memory=$1
	_threads=$2
	shift 2
	_holdids="$*"
	
	# Getting the options for the scheduler into variables
	# This is done so that the user can adapt these options to match them with his own scheduler

	#####################################################################################
	####### EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT #######
	#####################################################################################
	
	_OptExec=""					# Executable name for job scheduler
	_OptPriority=""					# Option for setting priority
	_OptMemory=""					# String for requesting memory
	_OptTerse=""					# Option for requesting output of a job id only
	_OptJoin=""					# String for joining the output files
	_OptOutput=""					# String for specifying output filename
	_OptHoldID=""					# String for requesting a job dependency
	_OptThreads=""					# String for requesting multiple processors on a node
	
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
		for _id in ${_holdids}
		do
			submitCMD="${submitCMD} ${_OptHoldID} ${_id}"
		done
	fi

	#####################################################################################
	####### STOP EDITING! STOP EDITING! STOP EDITING! STOP EDITING! STOP EDITING! #######
	#####################################################################################
	
	if [ -n "$submitCMD" ]
	then
		echo $submitCMD
	else
		echo -e "Job scheduler not defined! Please modify src/schedulerSettings/SettingsMISC.sh from the MUSE source package."
		exit 1
	fi
}

setDeleteJob()
{
	# Getting the job IDs
	_jobid="$*"
	
	# Getting the options for the scheduler into variables
	# This is done so that the user can adapt these options to match them with his own scheduler

	#####################################################################################
	####### EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT EDIT #######
	#####################################################################################

	_OptJobDeleteExec=""				# Executable name for deleting jobs
	
	#####################################################################################
	####### STOP EDITING! STOP EDITING! STOP EDITING! STOP EDITING! STOP EDITING! #######
	#####################################################################################

	echo ${_OptJobDeleteExec}
}
