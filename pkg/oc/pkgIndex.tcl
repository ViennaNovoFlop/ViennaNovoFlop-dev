# FILE: pkgIndex.tcl
#
# Last modified on: $Date: 2007-03-21 01:02:41 $
# Last modified by: $Author: donahue $
#
# Do not override an existing ifneeded script (from C, for example).
if {![string match "" [package ifneeded Oc 1.2.0.4]]} {return}
package ifneeded Oc 1.2.0.4 [list uplevel #0 [list source [file join $dir oc.tcl]]]
