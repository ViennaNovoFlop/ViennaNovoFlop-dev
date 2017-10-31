# unknown.tcl
#
# Configuration feature definitions for the configuration 'unknown'
#
# Editing instructions begin at "START EDIT HERE" below.

set config [Oc_Config RunPlatform]

set scriptfn [Oc_DirectPathname [info script]]
if {![string match [string tolower [file rootname [file tail $scriptfn]]] \
        [$config GetValue platform_name]]} {
    error "Configuration file '$scriptfn'
sourced by '[$config GetValue platform_name]'"
}

set localfn [file join [file dirname $scriptfn] local \
                [file tail $scriptfn]]
if {[file readable $localfn]} {
    if {[catch {source $localfn} msg]} {
        global errorInfo errorCode
	set msg [join [split $msg \n] \n\t]
	error "Error sourcing local platform file:\n    $localfn:\n\t$msg" \
		$errorInfo $errorCode
    }
}

########################################################################
# START EDIT HERE
# OOMMF software cannot classify your computing platform type.  See
# the Installation section of the OOMMF User Manual for instructions on 
# how to add a new platform type to the collection of types recognized 
# by OOMMF software.
#
# Say you add the new platform type 'foo' to describe your computing
# platform.  Then copy this file to ./foo.tcl , and edit it to
# describe your computing platform.  In its initial state, this file
# describes a rather generic Linux system, so there shouldn't be much
# editing required for any Unix-based platform.  Other platforms may
# be more difficult or impossible.  In any event, please e-mail the
# OOMMF developers for assistance setting up OOMMF for your particular
# circumstances.
########################################################################
# In order to properly build, install, and run on your computing 
# platform, the OOMMF software must know certain features of your
# computing environment.  In this file are lines which set the value of 
# certain features of your computing environment.  Each line looks like:
#
# $config SetValue <feature> {<value>}
#
# where each <feature> is the name of some feature of interest,
# and <value> is the value which is assigned to that feature in a
# description of your computing environment.  Your task is to edit
# the values as necessary to properly describe your computing 
# environment.
#
# The character '#' at the beginning of a line is a comment character.
# It causes the contents of that line to be ignored.  To select
# among lines providing alternative values for a feature, uncomment the
# line containing the proper value.
#
# The features in this file are divided into three sections.  The first
# section (REQUIRED CONFIGURATION) includes features which require you 
# to provide a value.  The second section (OPTIONAL CONFIGURATION)
# includes features which have usable default values, but which you
# may wish to customize.  The third section (ADVANCED CONFIGURATION)
# contains features which you probably do not need or want to change
# without a good reason.
########################################################################
# REQUIRED CONFIGURATION

########################################################################
# OPTIONAL CONFIGURATION

# Set the feature 'program_compiler_c++' to the program to run on this 
# platform to compile source code files written in the language C++ into 
# object files.  Select from the choices below.  If the compiler is not 
# in your path, be sure to use the whole pathname.  Also include any 
# options required to instruct your compiler to only compile, not link.  
#
# If your compiler is not listed below, additional features will
# have to be added in the ADVANCED CONFIGURATION section below to 
# describe to the OOMMF software how to operate your compiler.  Send
# e-mail to the OOMMF developers for assistance.
#
# The GNU C++ compiler 'g++'
# <URL:http://www.gnu.org/software/gcc/gcc.html>
# <URL:http://egcs.cygnus.com/>
$config SetValue program_compiler_c++ {g++ -c}

# Set the feature 'path_directory_temporary' to the name of an existing 
# directory on your computer in which OOMMF software should write 
# temporary files.  All OOMMF users must have write access to this 
# directory.
#
# $config SetValue path_directory_temporary {/tmp}

########################################################################
# ADVANCED CONFIGURATION

# Compiler option processing...
if {[string match g++ [file tail [lindex \
        [$config GetValue program_compiler_c++] 0]]]} {
    # ...for GNU g++ C++ compiler
    set opts [list -O%s]

    # Default warnings disable
    set nowarn [list -Wno-non-template-friend]
    if {[info exists nowarn] && [llength $nowarn]>0} {
       set opts [concat $opts $nowarn]
    }
    catch {unset nowarn}

    $config SetValue program_compiler_c++_option_opt "format \"$opts\""

    $config SetValue program_compiler_c++_option_out {format "-o \"%s\""}
    $config SetValue program_compiler_c++_option_src {format \"%s\"}
    $config SetValue program_compiler_c++_option_inc {format "\"-I%s\""}
    $config SetValue program_compiler_c++_option_def {format "\"-D%s\""}

    # Widest natively support floating point type
    $config SetValue program_compiler_c++_typedef_realwide "double"

    # Directories to exclude from explicit include search path, i.e.,
    # the -I list.  Some versions of gcc complain if "system" directories
    # appear in the -I list.
    $config SetValue \
    	program_compiler_c++_system_include_path [list /usr/include]
}

# The program to run on this platform to link together object files and
# library files to create an executable binary.
#
# Use the selected compiler to control the linking.
$config SetValue program_linker [lindex \
        [$config GetValue program_compiler_c++] 0]

# Linker option processing...
if {[string match g++ [file tail [lindex \
        [$config GetValue program_linker] 0]]]} {
    # ...for GNU g++ as linker
    $config SetValue program_linker_option_obj {format \"%s\"}
    $config SetValue program_linker_option_out {format "-o \"%s\""}
    $config SetValue program_linker_option_lib {format \"%s\"}
    $config SetValue program_linker_uses_-L-l {1}
}

# The program to run on this platform to create a single library file out
# of many object files.
$config SetValue program_libmaker {ar cr}

if {[string match ar [file tail [lindex \
        [$config GetValue program_libmaker] 0]]]} {
    # Option processing for ar
    $config SetValue program_libmaker_option_obj {format \"%s\"}
    $config SetValue program_libmaker_option_out {format \"%s\"}
}

# The absolute, native filename of the null device
$config SetValue path_device_null {/dev/null}

# A partial Tcl command (or script) which when completed by lappending
# a file name stem and evaluated returns the corresponding file name for
# an executable on this platform
$config SetValue script_filename_executable {format %s}

# A partial Tcl command (or script) which when completed by lappending
# a file name stem and evaluated returns the corresponding file name for
# an object file on this platform
$config SetValue script_filename_object {format %s.o}

# A partial Tcl command (or script) which when completed by lappending
# a file name stem and evaluated returns the corresponding file name for
# a static library on this platform
$config SetValue script_filename_static_library {format lib%s.a}

########################################################################
unset config
