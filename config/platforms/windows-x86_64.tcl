# FILE: wintel.tcl
#
# Configuration feature definitions for the configuration 'wintel'
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

if {[catch {$config GetValue program_compiler_c++_override}] \
       && ![catch {$config GetValue program_compiler_c++} _]} {
   # If program_compiler_c++ is set, but program_compiler_c++_override
   # is not, then assume user set the former instead of the latter,
   # and so copy the former to the latter to preserve the setting
   # across the setting of program_compiler_c++ in the "REQUIRED
   # CONFIGURATION" section below.
   $config SetValue program_compiler_c++_override $_
}

## Support for the automated buildtest scripts
if {[info exists env(OOMMF_BUILDTEST)] && $env(OOMMF_BUILDTEST)} {
   source [file join [file dirname [info script]] buildtest.tcl]
}


########################################################################
# START EDIT HERE
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

# NOTE: The rest of the REQUIRED CONFIGURATION is required only
# for building OOMMF software from source code.  If you downloaded
# a distribution with pre-compiled executables, no more configuration
# is required.
#
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
# Microsoft Visual C++
# <URL:http://msdn.microsoft.com/visualc/>
$config SetValue program_compiler_c++ {cl /c}
#

########################################################################
# OPTIONAL CONFIGURATION

# Set the feature 'path_directory_temporary' to the name of an existing
# directory on your computer in which OOMMF software should write
# temporary files.  All OOMMF users must have write access to this
# directory.
#
# $config SetValue path_directory_temporary {C:\temp}
# $config SetValue path_directory_temporary {C:\}

########################################################################
# SUPPORT PROCEDURES
#
# Load routines to guess the CPU, determine compiler version, and
# provide appropriate cpu-specific and compiler version-specific
# optimization flags.
source [file join [file dirname [Oc_DirectPathname [info script]]]  \
         cpuguess-wintel.tcl]

########################################################################
# LOCAL CONFIGURATION
#
# The following options may be defined in the platforms/local/wintel.tcl
# file:
#
## Specify whether or not to build in thread support.
## Thread support is included automatically if the tclsh interpreter used
## during the build process is threaded.  If you have a thread enabled
## tclsh, but don't want oommf_threads, override here.
# $config SetValue oommf_threads 0  ;# 1 to force threaded build,
#                                   ## 0 to force non-threaded build.
#
## Specify the number of default threads.  This is only meaningful
## for builds with thread support.
# $config SetValue thread_count 4  ;# Replace '4' with desired thread count.
#
## Use SSE intrinsics?  If so, specify level here.  Set to 0 to not use
## SSE intrinsics.  Leave unset to get the default (which may depend
## on the selected compiler).
# $config SetValue sse_level 2  ;# Replace '2' with desired level
#
## Override default C++ compiler.  Note the "_override" suffix
## on the value name.
# $config SetValue program_compiler_c++_override {icl /nologo /c /GX /GR}
#
## Processor architecture for compiling.  The default is "generic"
## which should produce an executable that runs on any cpu model for
## the given platform.  Optionally, one may specify "host", in which
## case the build scripts will try to automatically detect the
## processor type on the current system, and select compiler options
## specific to that processor model.  The resulting binary will
## generally not run on other architectures.
# $config SetValue program_compiler_c++_cpu_arch generic
#
## Variable type used for array indices, OC_INDEX.  This is a signed
## type which by default is sized to match the pointer width.  You can
## force the type by setting the following option.  The value should
## be a three item list, where the first item is the name of the
## desired (signed) type, the second item is the name of the
## corresponding unsigned type, and the third is the width of these
## types, in bytes.  It is assumed that both the signed and unsigned
## types are the same width, as otherwise significant code breakage is
## expected.  Example:
# $config SetValue program_compiler_c++_oc_index_type {__int64 {unsigned __int64} 8}
#
## For OC_INDEX type checks.  If set to 1, then various segments in
## the code are activated which will detect some array index type
## mismatches at link time.  These tests are not comprehensive, and
## will probably break most third party code, but may be useful during
## development testing.
# $config SetValue program_compiler_c++_oc_index_checks 1
#
###################
# Default handling of local defaults:
#
if {[catch {$config GetValue oommf_threads}]} {
   # Value not set in platforms/local/wintel.tcl,
   # so use Tcl setting.
   global tcl_platform
   if {[info exists tcl_platform(threaded)] \
          && $tcl_platform(threaded)} {
      $config SetValue oommf_threads 1  ;# Yes threads
   } else {
      $config SetValue oommf_threads 0  ;# No threads
   }
}
$config SetValue thread_count_auto_max 4 ;# Arbitrarily limit
## maximum number of "auto" threads to 4.
if {[catch {$config GetValue thread_count}]} {
   # Value not set in platforms/local/wintel.tcl, so try
   # to get value from environment:
   if {[info exists env(NUMBER_OF_PROCESSORS)]} {
      set processor_count $env(NUMBER_OF_PROCESSORS)
      set auto_max [$config GetValue thread_count_auto_max]
      if {$processor_count>$auto_max} {
         # Limit automatically set thread count to auto_max
         set processor_count $auto_max
      }
      $config SetValue thread_count $processor_count
   }
}

if {[catch {$config GetValue program_compiler_c++_override} compiler] == 0} {
    $config SetValue program_compiler_c++ $compiler
}

########################################################################
# ADVANCED CONFIGURATION

# Compiler option processing...
if {[catch {$config GetValue program_compiler_c++} ccbasename]} {
   set ccbasename {}  ;# C++ compiler not selected
} else {
   set ccbasename [file tail [lindex $ccbasename 0]]
}

# Microsoft Visual C++ compiler
if {[string match cl $ccbasename]} {
   set compilestr [$config GetValue program_compiler_c++]
   if {![info exists cl_version]} {
      set cl_version [GuessClVersion [lindex $compilestr 0]]
   }
   lappend compilestr /nologo /GR ;# /GR turns on RTTI
   set cl_major_version [lindex $cl_version 0]

   if {[lindex $cl_major_version 0]>7} {
      # The exception handling specification switch "/GX"
      # is deprecated in version 8.  /EHa enables C++
      # exceptions with SEH exceptions, /EHs enables C++
      # exceptions without SEH exceptions, and /EHc sets
      # extern "C" to default to nothrow.
      lappend compilestr /EHac
   } else {
      lappend compilestr /GX
   }
   $config SetValue program_compiler_c++ $compilestr
   unset compilestr

   # Optimization options for Microsoft Visual C++
   #
   # VC++ 6 (1998) and earlier are not supported.
   #
   # Options for VC++ 7.0 (2002), 7.1 (2003):
   #            Disable optimizations: /Od
   #             Maximum optimization: /Ox
   #      Enable runtime debug checks: /GZ
   #   Optimize for Pentium processor: /G5
   #         Optimize for Pentium Pro: /G6
   #
   # Options for VC++ 8.0 (2005), 9.0 (2008):
   #                  Disable optimizations: /Od
   #                   Maximum optimization: /Ox
   #                    Enable stack checks: /GZ
   #                   Require SSE2 support: /arch:SSE2
   # Fast (less predictable) floating point: /fp:fast
   #     Use portable but insecure lib fcns: /D_CRT_SECURE_NO_DEPRECATE
   #
   # Default optimization
   #   set opts {}
   # Max optimization
   set opts [GetClGeneralOptFlags $cl_version]
   # Aggressive optimization flags, some of which are specific to
   # particular cl versions, but are all processor agnostic.  CPU
   # specific opts are introduced in farther below.  See
   # cpuguess-wintel.tcl and x86-support.tcl for details.

   # CPU model architecture specific options.  To override, set value
   # program_compiler_c++_cpu_arch in
   # oommf/config/platform/local/wintel.tcl.
   if {[catch {$config GetValue program_compiler_c++_cpu_arch} cpu_arch]} {
      set cpu_arch generic
   }

   set cpuopts {}
   if {![string match generic [string tolower $cpu_arch]]} {
      # Arch specific build.  If cpu_arch is "host", then try to
      # guess.  Otherwise, assume user knows what he is doing and has
      # inserted an appropriate cpu_arch string, i.e., one that
      # matches the format and known types as returned from GuessCpu.
      if {[string match host $cpu_arch]} {
         set cpu_arch [GuessCpu]
      }
      # Use/don't use SSE intrinsics.  In the cpu_arch==host case,
      # the default behavior is to set this from the third element
      # of the GuessCpu return.  If cpu_arch!=host, then the
      # default is no.  You can always override the default
      # behavior setting the $config sse_level value in the local
      # platform file (see LOCAL CONFIGURATION above).
      if {[catch {$config GetValue sse_level}]} {
         # sse_level not set in LOCAL CONFIGURATION block
         $config SetValue sse_level [lindex $cpu_arch 2]
      }
      set cpuopts [GetClCpuOptFlags $cl_version $cpu_arch]
   }
   unset cpu_arch
   # You can override the above results by directly setting or
   # unsetting the cpuopts variable, e.g.,
   #
   #    set cpuopts [list /arch:SSE2]
   # or
   #    unset cpuopts
   #
   if {[info exists cpuopts] && [llength $cpuopts]>0} {
      set opts [concat $opts $cpuopts]
   }

   # Use/don't use SSE intrinsics.  The default is '2', because x86_64
   # guarantees at least SSE2.  You can override the value by setting
   # the $config sse_level value in the local platform file (see LOCAL
   # CONFIGURATION above).
   if {[catch {$config GetValue sse_level}]} {
      $config SetValue sse_level 2
   }

   # Silence security warnings
   if {$cl_major_version>7} {
      lappend opts /D_CRT_SECURE_NO_DEPRECATE
   }

   # NOTE: If you want good performance, be sure to edit ../options.tcl
   #  or ../local/options.tcl to include the line
   #    Oc_Option Add * Platform cflags {-def NDEBUG}
   #  so that the NDEBUG symbol is defined during compile.
   $config SetValue program_compiler_c++_option_opt "format \"$opts\""

   $config SetValue program_compiler_c++_option_out {format "\"/Fo%s\""}
   $config SetValue program_compiler_c++_option_src {format "\"/Tp%s\""}
   $config SetValue program_compiler_c++_option_inc {format "\"/I%s\""}
   $config SetValue program_compiler_c++_option_warn {
      format "/W4 /wd4505 /wd4702"
   }
   #   Warning C4505 is about removal of unreferenced local functions.
   # This seems to be a common occurrence when using templates with the
   # so-called "Borland" model.
   #   Warning C4702 is about unreachable code.  A lot of warnings of
   # this type are generated in the STL; I'm not even sure they are all
   # true.
   #
   # $config SetValue program_compiler_c++_option_debug {format "/MLd"}
   $config SetValue program_compiler_c++_option_debug {format "/Zi"}
   $config SetValue program_compiler_c++_option_def {format "\"/D%s\""}

   # Use OOMMF supplied erf() error function
   $config SetValue program_compiler_c++_property_no_erf 1

   # Use _hypot() in place of hypot()
   $config SetValue program_compiler_c++_property_use_underscore_hypot 1

   # Use _getpid() in place of getpid()
   $config SetValue program_compiler_c++_property_use_underscore_getpid 1

   # Widest natively support floating point type
   $config SetValue program_compiler_c++_typedef_realwide "double"

   # Visual C++ 6.0 does not support direct conversion from
   # unsigned __int64 to double.  If automatic detection doesn't
   # work, set cl_version directly to 5, 6, or 7, as appropriate.
   $config SetValue program_compiler_c++_uint64_to_double 0
   if {$cl_major_version>=7} {
      $config SetValue program_compiler_c++_uint64_to_double 1
   }

   # The program to run on this platform to create a single library file out
   # of many object files.
   # Microsoft Visual C++'s library maker
   $config SetValue program_libmaker {link /lib}
   # If your link doesn't accept the /lib option, try this instead:
   # $config SetValue program_libmaker {lib}
   $config SetValue program_libmaker_option_obj {format \"%s\"}
   $config SetValue program_libmaker_option_out {format "\"/OUT:%s\""}

   # The program to run on this platform to link together object files and
   # library files to create an executable binary.
   # Microsoft Visual C++'s linker
   $config SetValue program_linker {link}
   # $config SetValue program_linker {link /DEBUG} ;# For debugging
   $config SetValue program_linker_option_obj {format \"%s\"}
   $config SetValue program_linker_option_out {format "\"/OUT:%s\""}
   $config SetValue program_linker_option_lib {format \"%s\"}
   $config SetValue program_linker_option_sub {format "\"/SUBSYSTEM:%s\""}
   $config SetValue TCL_LIB_SPEC [$config GetValue TCL_VC_LIB_SPEC]
   $config SetValue TK_LIB_SPEC [$config GetValue TK_VC_LIB_SPEC]
   # Note: advapi32 is needed for GetUserName function in Nb package.
   $config SetValue TK_LIBS {user32.lib advapi32.lib}
   $config SetValue TCL_LIBS {user32.lib advapi32.lib}
   $config SetValue program_linker_uses_-L-l {0}
   $config SetValue program_linker_uses_-I-L-l {0}

   unset cl_version
   unset cl_major_version
}
catch {unset ccbasename}

# The absolute, native filename of the null device
$config SetValue path_device_null {nul:}

# A partial Tcl command (or script) which when completed by lappending
# a file name stem and evaluated returns the corresponding file name for
# an executable on this platform
$config SetValue script_filename_executable {format %s.exe}

# A partial Tcl command (or script) which when completed by lappending
# a file name stem and evaluated returns the corresponding file name for
# an object file on this platform
$config SetValue script_filename_object {format %s.obj}

# A partial Tcl command (or script) which when completed by lappending
# a file name stem and evaluated returns the corresponding file name for
# a static library on this platform
$config SetValue script_filename_static_library {format %s.lib}

# A list of partial Tcl commands (or scripts) which when completed by
# lappending a file name stem and evaluated returns the corresponding
# file name for an intermediate file produced by the linker on this platform
$config SetValue script_filename_intermediate [list \
   {format %s.ilk} {format %s.pdb} {format %s.map}]

########################################################################
unset config
