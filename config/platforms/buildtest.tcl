# File: buildtest.tcl
# Purpose: This script is used for automated building of
#  OOMMF with various build-time options, for regression
#  testing.
#
# The following environment variables *must* be set before
# sourcing this file:
#
#   OOMMF_BUILDTEST_COMPILER
#   OOMMF_BUILDTEST_CPUARCH
#   OOMMF_BUILDTEST_THREADS
#   OOMMF_BUILDTEST_NUMA
#
# Also, the "config" variable should be set to [Oc_Config RunPlatform]
#

# Append -warn 1 flag to cflags list, to enhance build diagnostics
if {[Oc_Option Get Platform cflags cf] == 0} {
   lappend cf -warn 1
} else {
   set cf [list -warn 1]
}
Oc_Option Add * Platform cflags $cf

$config SetValue program_compiler_c++_override $env(OOMMF_BUILDTEST_COMPILER)
$config SetValue program_compiler_c++_cpu_arch $env(OOMMF_BUILDTEST_CPUARCH)
$config SetValue oommf_threads $env(OOMMF_BUILDTEST_THREADS)
$config SetValue use_numa $env(OOMMF_BUILDTEST_NUMA)

# Signal that this file was loaded.  This is a safety check for runseq.tcl
set touchfile "[info script]-touch"
if {![file exists $touchfile]} {
   set chan [open $touchfile w]
   close $chan
}
if {[catch {file mtime $touchfile [clock seconds]}]} {
   # Prior to Tcl 8.3, 'file mtime' did not support changing the mtime.
   # This is also protection against 'file exists' being broken on some
   # file systems.
   set chan [open $touchfile w]
   close $chan
}
