# FILE: pkgIndex.tcl
#
# Last modified on: $Date: 2007-03-21 01:17:05 $
# Last modified by: $Author: donahue $
#
package ifneeded Net 1.2.0.4 [list uplevel #0 [list source [file join $dir net.tcl]]]
