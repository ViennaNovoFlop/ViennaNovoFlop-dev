#!/bin/sh
# FILE: cpuguess-lintel.tcl
#
# Wrapper Tcl script to guess the CPU model.
#
#    v--Edit here if necessary \
exec tclsh "$0" ${1+"$@"}
########################################################################

proc Usage {} {
   puts stderr "Usage: tclsh cpuguess.tcl <platform> ?extras?"
   puts stderr " where <platform> is, for example, \"lintel\","
   puts stderr " and optional parameter extras is 0 or 1."
   exit 1
}

if {[llength $argv] <1 || [llength $argv]>2} { Usage }
if {[lsearch -regexp $argv {^-*[h|H]}]>=0} { Usage }

set platform_file "cpuguess-"
append platform_file [file rootname [file tail [lindex $argv 0]]]
append platform_file ".tcl"
set platform_file [file join [file dirname [info script]]  $platform_file]

set extras 0
if {[llength $argv]>1} { set extras [lindex $argv 1] }

if {![file exists $platform_file]} {
   # puts stderr "Platform cpu file \"$platform_file\" does not exist"
   # exit 2
   puts "CPU guess: unknown"
   exit
}

if {![file readable $platform_file]} {
   puts stderr \
      "Platform cpu file \"$platform_file\" exists but is not readable"
   exit 2
}

source $platform_file

set cpu_arch [GuessCpu]
puts "CPU guess: $cpu_arch"

if {$extras} {
   foreach elt [info proc Guess*Version] {
      if {[regexp -- {^Guess(.*)Version$} $elt dum compiler_str]} {
         set compiler [string tolower $compiler_str]
         if {![catch {$elt $compiler} compiler_version] \
                && ![string match {} $compiler_version]} {
            puts "$compiler version: $compiler_version"
            set genopts [Get${compiler_str}GeneralOptFlags $compiler_version]
            set cpuopts [Get${compiler_str}CpuOptFlags \
                            $compiler_version $cpu_arch]
            puts "  $compiler flags: [concat $genopts $cpuopts]"
         }
      }
   }
}


