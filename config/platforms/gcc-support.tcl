# FILE: gcc-support.tcl
#
# A couple of Tcl procs for use determining compiler options for gcc:
#    GuessGccVersion
# which runs gcc to guess the version of GCC being used.
#    GetGccGeneralOptFlags
# which returns a list of aggressive, non-processor-specific
# optimization flags that can be passed to gcc.
#

# Routine to guess the gcc version.  The import, gcc, is used via "exec
# $gcc --version" (or, rather, the "open" analogue) to determine the
# gcc version (since the flags accepted by gcc vary by version).
# Return value is the gcc version string as a list of numbers,
# for example, {4 1 1} for version 4.1.1.
proc GuessGccVersion { gcc } {
    set guess {}
    catch {
	set fptr [open "|$gcc --version" r]
	set verstr [read $fptr]
	close $fptr
        set digstr {[0-9]+\.[0-9.]+}
        set ws "\[ \t\n\r\]"
	regexp -- "(^|$ws)($digstr)($ws|$)" $verstr dummy dum0 guess dum1
    }
    return [split $guess "."]
}

# Routine that returns aggressive, processor agnostic, optimization
# flags for gcc.  The import is the gcc version, as returned by the
# preceding GuessGccVersion proc.  Note that the flags accepted by gcc
# vary by version.
#
proc GetGccGeneralOptFlags { gcc_version } {

   if {[llength $gcc_version]<2} {
      # Unable to determine gcc version.  Return an empty string.
      return {}
   }
   set verA [lindex $gcc_version 0]
   set verB [lindex $gcc_version 1]

   if {![regexp -- {[0-9]+} $verA] || ![regexp -- {[0-9]+} $verB]} {
      return -code error "Invalid input:\
         gcc_version should be a list of numbers, not\
         \"$gcc_version\""
   }

   # Optimization options
   set opts [list -O2 -ffast-math]
   if {$verA>=3} {
      lappend opts -frename-registers

      # Docs say -fstrict-aliasing is in 2.95.  I can't find docs
      # for 2.8, but one extant installation I have (on SGI) doesn't
      # recognize -fscrict-aliasing.  To be safe, just require 3.x.
      # OTOH, i686-apple-darwin9-gcc-4.0.1 (OS X Leopard) (others?)
      # has an optimization bug yielding a "Bus error" for some
      # code if -fstrict-aliasing is enabled.
      if {$verA!=4 || $verB>0} {
         lappend opts -fstrict-aliasing
      }

      if {$verA>=4 || $verB>=4} {
         lappend opts -fweb
      }
   }

   # Frame pointer: Some versions of gcc don't handle exceptions
   # properly w/o frame-pointers.  This typically manifests as
   # Oxs dying without an error message while loading a MIF file.
   # On some x86 systems, including in addition the flag
   # -momit-leaf-frame-pointer makes the problem go away.
   # (See proc GetGccCpuOptFlags in x86-support.tcl).
   # Comment this out if the aforementioned problem occurs.
   lappend opts -fomit-frame-pointer

   return $opts
}
