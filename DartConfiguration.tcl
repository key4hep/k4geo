# This file is configured by CMake automatically as DartConfiguration.tcl
# If you choose not to use CMake, this file may be hand configured, by
# filling in the required variables.


# Configuration directories and files
SourceDirectory: /project/lepcol/users/cligtenb/ilcsoftPixel/lcgeo
BuildDirectory: /project/lepcol/users/cligtenb/ilcsoftPixel/lcgeo

# Where to place the cost data store
CostDataFile: 

# Site is something like machine.domain, i.e. pragmatic.crd
Site: drac.nikhef.nl

# Build name is osname-revision-compiler, i.e. Linux-2.4.2-2smp-c++
BuildName: Linux-2.6.32-696.1.1.el6.x86_64-/project/lepcol/ilcsoft/dependencies/gcc-4.9.4/bin/g++-RelWithDebInfo

# Submission information
IsCDash: TRUE
CDashVersion: 
QueryCDashVersion: 
DropSite: aidasoft.desy.de
DropLocation: /CDash/submit.php?project=lcgeo
DropSiteUser: 
DropSitePassword: 
DropSiteMode: 
DropMethod: http
TriggerSite: 
ScpCommand: /usr/bin/scp

# Dashboard start time
NightlyStartTime: 01:00:00 UTC

# Commands for the build/test/submit cycle
ConfigureCommand: "/project/lepcol/ilcsoft/dependencies/CMake/3.6.3/bin/cmake" "/project/lepcol/users/cligtenb/ilcsoftPixel/lcgeo"
MakeCommand: /project/lepcol/ilcsoft/dependencies/CMake/3.6.3/bin/cmake --build . --config "${CTEST_CONFIGURATION_TYPE}"
DefaultCTestConfigurationType: Release

# version control
UpdateVersionOnly: 

# CVS options
# Default is "-d -P -A"
CVSCommand: /usr/bin/cvs
CVSUpdateOptions: -d -A -P

# Subversion options
SVNCommand: /usr/bin/svn
SVNOptions: 
SVNUpdateOptions: 

# Git options
GITCommand: /project/lepcol/ilcsoft/dependencies/git2/bin/git
GITInitSubmodules: 
GITUpdateOptions: 
GITUpdateCustom: 

# Perforce options
P4Command: P4COMMAND-NOTFOUND
P4Client: 
P4Options: 
P4UpdateOptions: 
P4UpdateCustom: 

# Generic update command
UpdateCommand: /project/lepcol/ilcsoft/dependencies/git2/bin/git
UpdateOptions: 
UpdateType: git

# Compiler info
Compiler: /project/lepcol/ilcsoft/dependencies/gcc-4.9.4/bin/g++
CompilerVersion: 4.9.4

# Dynamic analysis (MemCheck)
PurifyCommand: 
ValgrindCommand: 
ValgrindCommandOptions: 
MemoryCheckType: 
MemoryCheckSanitizerOptions: 
MemoryCheckCommand: /usr/bin/valgrind
MemoryCheckCommandOptions: 
MemoryCheckSuppressionFile: 

# Coverage
CoverageCommand: /project/lepcol/ilcsoft/dependencies/gcc-4.9.4/bin/gcov
CoverageExtraFlags: -l

# Cluster commands
SlurmBatchCommand: SLURM_SBATCH_COMMAND-NOTFOUND
SlurmRunCommand: SLURM_SRUN_COMMAND-NOTFOUND

# Testing options
# TimeOut is the amount of time in seconds to wait for processes
# to complete during testing.  After TimeOut seconds, the
# process will be summarily terminated.
# Currently set to 25 minutes
TimeOut: 1500

# During parallel testing CTest will not start a new test if doing
# so would cause the system load to exceed this value.
TestLoad: 

UseLaunchers: 
CurlOptions: 
# warning, if you add new options here that have to do with submit,
# you have to update cmCTestSubmitCommand.cxx

# For CTest submissions that timeout, these options
# specify behavior for retrying the submission
CTestSubmitRetryDelay: 5
CTestSubmitRetryCount: 3
