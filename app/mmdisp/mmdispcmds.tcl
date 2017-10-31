# FILE: mmdispcmds.tcl
#
#	Pure Tcl commands supporting mmDisp
#
# Last modified on: $Date: 2004-09-02 13:01:41 $
# Last modified by: $Author: dgp $

# Verify that C++ portion of this version of the Mmdispcmds extension 
# has been initialized
#
# NOTE: version number below must match that in ./mmdispcmds.h
package require -exact Mmdispcmds 1.2.0.4

Oc_CheckTclIndex Mmdispcmds

# Set up for autoloading of Mmdispcmds commands implemented in Tcl
set mmdispcmds(library) [file dirname [info script]]
if { [lsearch -exact $auto_path $mmdispcmds(library)] == -1 } {
    lappend auto_path $mmdispcmds(library)
}

