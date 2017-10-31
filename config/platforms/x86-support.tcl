# FILE: x86-support.tcl
#
# Several Tcl procs for use on x86-based systems:
#    GuessCpuArch_VendorFamilyModel
# which returns a cpu-arch list, {vendor type sse-level} given the
# vendor, processor family and model id numbers.
#    GuessCpuArch_NameStr
# which produces a cpu-arch list like the previous proc, but using
# the "model name" string.  (Avoid using this proc.)
#    GetGccCpuOptFlags
# which returns a list of processor-specific optimization flags that
# can be passed to gcc.
#
# Additionally, this file sources gcc-support.tcl, which provides
#    GuessGccVersion
# which runs gcc to guess the version of GCC being used.
#    GetGccGeneralOptFlags
# which returns a list of aggressive, non-processor-specific
# optimization flags that can be passed to gcc.
#
# NOTE: I don't know any definitive way to ascertain sse level for
#  the AMD k8 chips (i.e., sse2 versus sse3) based on just the
#  family and model numbers.  On linux-based systems, it is safer
#  and easier to get this information from the "flags" output of
#  /proc/cpuinfo ("pni" (Prescott New Instructions) means sse3) than
#  the cpu-arch result from either GuessCpuArch routine.
#
########################################################################

########################################################################
### Generic gcc support
source [file join [file dirname [info script]] gcc-support.tcl]


########################################################################
# Routines to guess the CPU architecture. Output is a three element list.
# The first element is the vendor, which will be one of:
#     unknown, intel, or amd.
# The second element is the cpu type, which will be one of:
#     unknown,
#     i386, i486, pentium, pentium-mmx, pentiumpro, pentium2,
#     pentium3, pentium-m, core, core2
#     pentium4, prescott, nocona, pentium-d,
#     k5, k6, k6-2, k6-3, athlon, athlon-4, k8
# The third element is the supported SSE level, which will be one of:
#     0, 1, 2, 3
#

proc GuessCpuArch_VendorFamilyModel { vendor cpu_family cpu_model } {
   set vendor  [string tolower $vendor]
   if {[string match "*intel*" $vendor]} {
      set vendor intel
   } elseif {[string match "*amd*" $vendor]} {
      set vendor amd
   } else {
      set vendor unknown
   }

   set cputype unknown
   set sselevel 0

   if {[string match intel $vendor]} {
      if {$cpu_family < 4} {
         set cputype i386
      } elseif {$cpu_family == 4} {
         set cputype i486
      } elseif {$cpu_family == 5} {
         set cputype pentium
         if {$cpu_model>3} { set cputype pentium-mmx }
      } elseif {$cpu_family == 6} {
         if {$cpu_model < 3} {
            set cputype pentiumpro
         } elseif {$cpu_model < 7} {
            set cputype pentium2
         } elseif {$cpu_model <9 \
                      || $cpu_model == 10 || $cpu_model == 11 } {
            set cputype pentium3
            set sselevel 1
         } elseif {$cpu_model == 9 || $cpu_model == 13} {
            set cputype pentium-m
            set sselevel 2
         } elseif {$cpu_model == 14} {
            set cputype core
            set sselevel 3
         } elseif {$cpu_model >= 15} {
            set cputype core2
            set sselevel 3
         }
      } elseif {$cpu_family >= 15} {
         set cputype pentium4
         set sselevel 2
         if {$cpu_model >=3} {
            set sselevel 3
            set cputype prescott
            if {$cpu_model >= 4} { set cputype nocona    }
            if {$cpu_model >= 6} { set cputype pentium-d }
         }
      }
   } elseif {[string match amd $vendor]} {
      set vendor "amd"
      if {$cpu_family == 4} {
         if {$cpu_model < 14} {
            set cputype i486
         } else {
            set cputype pentium
         }
      } elseif {$cpu_family == 5} {
         if {$cpu_model < 6} {
            set cputype k5
         } elseif {$cpu_model < 8} {
            set cputype k6
         } elseif {$cpu_model == 8} {
            set cputype k6-2
         } else {
            set cputype k6-3
         }
      } elseif {$cpu_family == 6} {
         if {$cpu_model<6} {
            set cputype athlon
         } else {
            set cputype athlon-4
            set sselevel 1
         }
      } elseif {$cpu_family >= 15} {
         set cputype k8
         set sselevel 2
         if {$cpu_model > 32 } {  ;# Best guess?
            set sselevel 3
         }
      }
   }
   return [list $vendor $cputype $sselevel]
}

proc GuessCpuArch_NameStr { namestr } {
   # Seriously, don't use this routine.  It is inferior in every way
   # to GuessCpuArch_VendorFamilyModel.  Some common names like "Xeon"
   # and "Sempron" are entirely missing, because on their own "Xeon"
   # and "Sempron" carry very little meaning.
   #   This proc is hanging around on the off chance that someone
   # might find some use for it someday.

   set namestr [string tolower $namestr]
   if {[string match "*intel*" $namestr]} {
      set vendor intel
   } elseif {[string match "*amd*" $namestr]} {
      set vendor amd
   } else {
      set vendor unknown
   }

   set sselevel 0
   set cputype unknown
   switch -regexp -- $namestr {
      {^386}                    { set cputype i386 }
      {^486}                    { set cputype i486 }
      {intel(\(r\)|) *core(\(tm\)|) *2}  {
         set cputype core2
         set sselevel 3
      }
      {intel(\(r\)|) *core}     {
         set cputype core
         set sselevel 3
      }
      {pentium(\(r\)|) *d( |$)}  {
         set cputype pentium-d
         set sselevel 3
      }
      {pentium(\(r\)|) *(iv|4)}  {
         set cputype pentium4
         set sselevel 2
      }
      {pentium(\(r\)|) *m( |$)}  {
         set cputype pentium-m
         set sselevel 2
      }
      {pentium(\(r\)|) *(iii|3)} {
         set cputype pentium3
         set sselevel 1
      }
      {pentium(\(r\)|) *(ii|2)} { set cputype pentium2 }
      {pentium(\(r\)|) *pro}    { set cputype pentiumpro }
      {pentium}                 { set cputype pentium }
      {overdrive podp5v83}      { set cputype pentium }
      {opteron}                 -
      {athlon(\(tm\)|) 64}      -
      {athlon(\(tm\)|) fx}      {
         set cputype k8
         set sselevel 2
      }
      {athlon(\(tm\)|) mp}      -
      {athlon(\(tm\)|) xp}      -
      {athlon(\(tm\)|) (iv|4)}  {
         set cputype athlon-4
         set sselevel 1
      }
      {athlon(\(tm\)|) (tbird)} -
      {athlon}                  { set cputype athlon }
      {^k8}                     {
         set cputype k8
         set sselevel 2
      }
      {^k6-3}                   { set cputype k6-3 }
      {^k6-2}                   { set cputype k6-2 }
      {^k6}                     { set cputype k6 }
      {^k5}                     { set cputype k5 }
   }

   return [list $vendor $cputype $sselevel]
}

# Routine that determines processor specific optimization flags for
# gcc.  The first import is the gcc version, as returned by the
# GuessGccVersion proc (see gcc-support.tcl).  Note that the flags
# accepted by gcc vary by version.  The second import, cpu_arch,
# should match output from the GuessCpu proc above.  Return value is a
# list of gcc flags.
proc GetGccCpuOptFlags { gcc_version cpu_arch } {
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

   # Extract cpu information from cpu_arch import
   set cpu_vendor "unknown"
   set cpu_type   "unknown"
   set cpu_sse    0
   foreach {cpu_vendor cpu_type cpu_sse} $cpu_arch { break }

   # Construct optimization flags
   # Note: -fprefetch-loop-arrays is available in gcc 3.1
   # and later, for Intel processors pentium3 or better,
   # and for AMD processors k6-2 or better.
   set cpuopts {}
   if {$verA==2 && $verB>=95} {
      # Don't bother setting -march in case of
      # i386, i486, or k5
      switch -glob -- $cpu_type {
         k6*         { set cpuopts [list -march=k6] }
         pentium  -
         pentium-mmx { set cpuopts [list -march=pentium] }
         pentium* -
         prescott -
         nocona   -
         core*    -
         athlon*  -
         opteron  -
         k8          { set cpuopts [list -march=pentiumpro] }
      }
      set cpu_sse 0
   } elseif {$verA==3 && $verB<1} {
      # Don't bother setting -march in case of
      # i386, i486, or k5
      switch -glob -- $cpu_type {
         pentium  -
         pentium-mmx { set cpuopts [list -march=pentium] }
         pentium* -
         prescott -
         nocona   -
         core*       { set cpuopts [list -march=pentiumpro] }
         k6*         { set cpuopts [list -march=k6] }
         athlon*  -
         opteron  -
         k8          { set cpuopts [list -march=athlon] }
      }
      set cpu_sse 0
   } elseif {$verA==3  && $verB<3} {
      # Don't bother setting -march in case of
      # i386, i486, or k5
      switch -glob -- $cpu_type {
         pentium    -
         pentium-mmx { set cpuopts [list -march=pentium] }
         pentiumpro { set cpuopts [list -march=pentiumpro] }
         pentium2   { set cpuopts [list -march=pentium2] }
         pentium3   -
         pentium-m  { set cpuopts [list -march=pentium3 \
                                      -fprefetch-loop-arrays] }
         pentium*   -
         prescott   -
         nocona     -
         core*      { set cpuopts [list -march=pentium4 \
                                      -fprefetch-loop-arrays] }
         k6         { set cpuopts [list -march=k6] }
         k6-2       { set cpuopts [list -march=k6-2 -fprefetch-loop-arrays] }
         k6-3       { set cpuopts [list -march=k6-3 -fprefetch-loop-arrays] }
         athlon     { set cpuopts [list -march=athlon -fprefetch-loop-arrays] }
         athlon-tbird { set cpuopts [list -march=athlon-tbird \
                                        -fprefetch-loop-arrays] }
         athlon*    -
         opteron    -
         k8         { set cpuopts [list -march=athlon-4 \
                                      -fprefetch-loop-arrays] }
      }
      if {$cpu_sse>=2} { set cpu_sse 2 }
   } elseif {$verA==3 && $verB==3} {
      # Don't bother setting -march in case of
      # i386, i486, or k5
      switch -glob -- $cpu_type {
         pentium     -
         pentium-mmx { set cpuopts [list -march=pentium] }
         pentiumpro  { set cpuopts [list -march=pentiumpro] }
         pentium2    { set cpuopts [list -march=pentium2] }
         pentium3    -
         pentium-m   { set cpuopts [list -march=pentium3 \
                                       -fprefetch-loop-arrays] }
         pentium4    { set cpuopts [list -march=pentium4 \
                                       -fprefetch-loop-arrays] }
         prescott    { set cpuopts [list -march=prescott \
                                       -fprefetch-loop-arrays] }
         nocona      -
         pentium*    -
         core*       { set cpuopts [list -march=nocona \
                                       -fprefetch-loop-arrays] }
         k6          { set cpuopts [list -march=k6] }
         k6-2        { set cpuopts [list -march=k6-2 -fprefetch-loop-arrays] }
         k6-3        { set cpuopts [list -march=k6-3 -fprefetch-loop-arrays] }
         athlon      { set cpuopts [list -march=athlon \
                                       -fprefetch-loop-arrays] }
         athlon-tbird { set cpuopts [list -march=athlon-tbird \
                                        -fprefetch-loop-arrays] }
         athlon*     -
         opteron     -
         k8          { set cpuopts [list -march=athlon-4 \
                                       -fprefetch-loop-arrays] }
      }
      if {$cpu_sse>=3} { set cpu_sse 3 }  ;# Safety
   } elseif {($verA==3 && $verB>=4) || ($verA==4 && $verB<=1)} {
      # Don't bother setting -march in case of
      # i386, i486, or k5
      switch -glob -- $cpu_type {
         pentium     -
         pentium-mmx { set cpuopts [list -march=pentium] }
         pentiumpro  { set cpuopts [list -march=pentiumpro] }
         pentium2    { set cpuopts [list -march=pentium2] }
         pentium3    { set cpuopts [list -march=pentium3 \
                                       -fprefetch-loop-arrays] }
         pentium-m   { set cpuopts [list -march=pentium-m \
                                       -fprefetch-loop-arrays] }
         pentium4    { set cpuopts [list -march=pentium4 \
                                       -fprefetch-loop-arrays] }
         prescott    { set cpuopts [list -march=prescott \
                                       -fprefetch-loop-arrays] }
         nocona      -
         pentium*    -
         core*       { set cpuopts [list -march=nocona \
                                       -fprefetch-loop-arrays] }
         k6          { set cpuopts [list -march=k6 -fprefetch-loop-arrays] }
         k6-2        { set cpuopts [list -march=k6-2 -fprefetch-loop-arrays] }
         k6-3        { set cpuopts [list -march=k6-3 -fprefetch-loop-arrays] }
         athlon      { set cpuopts [list -march=athlon \
                                       -fprefetch-loop-arrays] }
         athlon-4    { set cpuopts [list -march=athlon-4 \
                                       -fprefetch-loop-arrays] }
         opteron     -
         k8          { set cpuopts [list -march=k8 -fprefetch-loop-arrays] }
      }
      if {$cpu_sse>=3} { set cpu_sse 3 }  ;# Safety
   } elseif {$verA>4 || ($verA==4 && $verB>=2)} {
      set cpuopts [list -march=native]
      # -march/-mtune=native setting introduced with gcc 4.2
      switch -glob -- $cpu_type {
         pentium     -
         pentium-mmx -
         pentiumpro  -
         pentium2    {}
         pentium*    -
         prescott    -
         nocona      -
         pentium*    -
         core*       { lappend cpuopts -fprefetch-loop-arrays }
         k6          {}
         k6-*        -
         athlon*     -
         opteron     -
         k8          { lappend cpuopts -fprefetch-loop-arrays }
      }
      if {$cpu_sse>=3} { set cpu_sse 3 }  ;# Safety
   }
   if {$cpu_sse>0} {
      lappend cpuopts -mfpmath=sse -msse
      for {set sl 2} {$sl<=$cpu_sse} {incr sl} {
         lappend cpuopts -msse$sl
      }
   }

   # Frame pointer: Some versions of gcc don't handle exceptions
   # properly w/o frame-pointers.  This typically manifests as
   # Oxs dying without an error message while loading a MIF file.
   # Interestingly, including -momit-leaf-frame-pointer appears
   # to work around this problem, at least on some systems.  YMMV;
   # Comment this out if the aforementioned problem occurs.
   lappend cpuopts -momit-leaf-frame-pointer

   return $cpuopts
}

# Routine to guess the Intel C++ version.  The import, icpc, is used
# via "exec $icpc --version" (or, rather, the "open" analogue) to
# determine the icpc version (since the flags accepted by icpc vary by
# version).  Return value is the icpc version string as a list of
# numbers, for example, {10 1} for version 10.1
proc GuessIcpcVersion { icpc } {
    set guess {}
    catch {
	set fptr [open "|$icpc --version" r]
	set verstr [read $fptr]
	close $fptr
        set digstr {[0-9]+\.[0-9.]+}
        set ws "\[ \t\n\r\]"
	regexp -- "(^|$ws)($digstr)($ws|$)" $verstr dummy dum0 guess dum1
    }
    return [split $guess "."]
}

# Routines that report optimization flags for icpc.  Right now these
# are mostly placeholders, but they may be expanded and refined in the
# future.
proc GetIcpcGeneralOptFlags { icpc_version } {
   set opts [list -O3 -ipo -no-prec-div -ansi_alias \
                -fp-model fast=2 -fp-speculation fast]
   return $opts
}
proc GetIcpcCpuOptFlags { icpc_version cpu_arch } {
   # CPU model architecture specific options.
   set cpuopts {}
   set cpu_vendor "unknown"
   set cpu_type   "unknown"
   set cpu_sse    0
   foreach {cpu_vendor cpu_type cpu_sse} $cpu_arch {
      set cpu_vendor [string tolower $cpu_vendor]
      set cpu_type   [string tolower $cpu_type]
      break
   }
   # The following flags are based on the docs for the icc version
   # 10.0.  Allowed flags may differ for other versions of icc.
   # Note: -xT and -xP imply SSE3.
   switch -exact -- $cpu_type {
      athlon-4  { lappend cpuopts -xK }
      opteron   -
      k8        { lappend cpuopts -xW}
      pentium3  { lappend cpuopts -xK }
      pentium-m { lappend cpuopts -xW}
      pentium4  { lappend cpuopts -xN}
      prescott  -
      nocona    -
      pentium-d -
      core      { lappend cpuopts -xP }
      core2     { lappend cpuopts -xT }
      default   {
         if {$cpu_sse>=3 && [string match intel $cpu_vendor]} {
            lappend cpuopts -xP
         } elseif {$cpu_sse>=2} {
            lappend cpuopts -xW
         } elseif {$cpu_sse>=1} {
            lappend cpuopts -xK
         }
      }
   }

   return $cpuopts
}
