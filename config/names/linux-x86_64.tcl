# linux-x86_64.tcl
#
# Defines the Oc_Config name 'linux-x86_64' to indicate the Linux
# operating system running on the x86_64 architecture.

Oc_Config New _ [string tolower [file rootname [file tail [info script]]]] {
    global tcl_platform env
    if {![regexp -nocase -- linux $tcl_platform(os)]} {
        return 0
    }
    if {![string match x86_64 $tcl_platform(machine)]} {
	return 0
    }
    return 1
}
