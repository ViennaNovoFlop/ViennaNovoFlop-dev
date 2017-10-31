# wintel.tcl
#
# Defines the Oc_Config name 'wintel' to indicate the Windows 95/NT operating
# system running on an Intel architecture.

Oc_Config New _ [string tolower [file rootname [file tail [info script]]]] {
   global tcl_platform
   if {![regexp -nocase -- windows $tcl_platform(platform)]} {
      return 0   ;# Not windows
   }
   if {([info exists tcl_platform(pointerSize)] &&
        $tcl_platform(pointerSize) == 8) ||
       (![info exists tcl_platform(pointerSize)] &&
        [string match amd64 $tcl_platform(machine)])} {
      # Looks like 64-bit Windows.  Just check that this
      # isn't a cygwin variant.
      return [expr {![Oc_IsCygwinPlatform]}]
   }
   return 0   ;# Not x86_64
}
