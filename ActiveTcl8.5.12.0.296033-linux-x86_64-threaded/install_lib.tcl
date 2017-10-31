# -*- tcl -*-
# Library to support the ActiveTcl installer scripts.
# Contains the common functionality. Its first action is to load the
# data file containing distribution/release dependent information.

package require Tcl 8.4

source [file join $::SCRIPT_DIR install_data.tcl]

set ::AT(MVERSION) $::AT(VERSION)
if {$::tcl_platform(machine) eq "amd64"} {
    # In 8.5 we could use pointerSize == 8
    # Ideally this would check what variant is being installed
    append ::AT(MVERSION) " (64-bit)"
}
if {$::AT(MATURITY) ne "stable" && $::AT(MATURITY) ne "final"} {
    append ::AT(MVERSION) " $::AT(MATURITY)"
}

# Compatibility code to make this library useable for already
# distributed releases.

if {![info exists ::AT(tcl_VERSION)]} {set ::AT(tcl_VERSION) [info tclversion]}
if {![info exists ::AT(DEBUG)]}       {set ::AT(DEBUG) 0}
if {![info exists ::AT(TDK_ATVERS)]}  {set ::AT(TDK_ATVERS) 8.3.5.0}

set AT(ND_VER) [string map {. ""} $::AT(tcl_VERSION)]
set AT(ND_VERG) $AT(ND_VER)
if {$AT(DEBUG)} {
    append AT(ND_VERG) g
    append AT(dPGROUP) " Debug"
}
set wish_exe wish$::AT(ND_VERG).exe


## Note aside: On unix we are patching the starkits to contain the
## absolute path of the installed tclsh. This means that the tools
## will work even if the PATH is not set. On windows we simply rely on
## the association for '.tcl'.

set INSTALL_MODE admin ; # Common installation is default.

# ======================================================================

proc pathext_possible {} {
    set ap [check_pathext_possible]
    proc pathext_possible {} [list return $ap]
    return $ap
}
proc check_pathext_possible {} {
    # We try to find out if parts of "admin" installation is possible.
    # We do this by execution a number of write-operations to the
    # registry which should fail if we don't have the required privileges.

    if {[string equal "Windows 95" $::tcl_platform(os)]} {
	return 1
    }

    global systemEnv

    if {[catch {set current [r_get $systemEnv PATHEXT]}]} {
	# Unable to read system information we need
	return 0
    } elseif {[catch {r_set $systemEnv PATHEXT "$current" sz}]} {
	# Unable to write into system information. (We are writing the
	# previously read value, this makes sure that we don't destroy
	# the data during this test).
	return 0
    }

    return 1
}

# Self-modifying code, caches result of complex check in the proc
# definition.

proc admin_possible {} {
    set ap [check_admin_possible]
    proc admin_possible {} [list return $ap]
    return $ap
}
proc check_admin_possible {} {
    # We try to find out if "admin" installation is possible.
    # We do this by execution a number of write-operations to the
    # registry which should fail if we don't have the required privileges.

    if {[string equal "Windows 95" $::tcl_platform(os)]} {
	return 1
    }

    global systemEnv

    if {[catch {set current [r_get $systemEnv PATH]}]} {
	# Unable to read system information we need
	return 0
    } elseif {[catch {r_set $systemEnv PATH "$current" expand_sz}]} {
	# Unable to write into system information. (We are writing the
	# previously read value, this makes sure that we don't destroy
	# the data during this test).
	return 0
    }

    return 1
}

global systemEnv ; set systemEnv {HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment}
global userEnv   ; set userEnv   {HKEY_CURRENT_USER\Environment}

# ======================================================================

# Return license text
#

# FUTURE -- Have the packaging code provide the installer with the
# name of README and License files, instead of having it to search for
# them. Via install_data.

proc license_file {} {
    global LICENSE AT

    set lf [file join $::SCRIPT_DIR $AT(LICENSE)]
    if {![file exists $lf]} {
	# Broken distribution.
	return -code error "Distribution broken.\n\
		License file \"$AT(LICENSE)\" not found."
    }
    set LICENSE [file tail $lf]
    return $lf
}

proc license_text {} {
    global LICENSE ; # For use in patch_registry.

    set fh [open [license_file] r]
    set data [read $fh]
    close $fh
    return $data
}

# Return overinstall warning
#
proc overinstall_warning {} {
    return "WARNING:
You are about to install $::AT(dPGROUP) $::AT(MVERSION) in a directory
containing an installation of $::AT(dPGROUP) v$::AT(InstVersion).

Choose 'Next' to install over the existing files.
Choose 'Back' to select a different installation directory."
}

# Return message on finish
#
# unix		\/ normal, lite, pro
# windows	/\

array set ::PARTING {}

set ::PARTING(unix) {Please do not forget to extend your PATH and MANPATH variables to
get access to the applications and manpages distributed with @name@.

For a csh or compatible perform
    setenv PATH "@idir@/bin:$PATH"

For a sh or similar perform
    PATH="@idir@/bin:$PATH"
    export PATH

Some shells (bash for example) allow
    export PATH="@idir@/bin:$PATH"

Similar changes are required for MANPATH}

set ::PARTING(windows) {}

proc parting_message {} {
    # The parting message depends on the platform.

    set key $::tcl_platform(platform)

    if {![info exists ::RUNTIME_DIR]} {
	set ::RUNTIME_DIR "<installDirectory>"
    }
    return [string map [list \
	    @name@ $::AT(NAME) \
	    @idir@ $::RUNTIME_DIR \
	    ] $::PARTING($key)]
}

proc hype_message {} {
    if {![string equal $::AT(MODE) normal]} {return}
    log  {  }                                                                hype
    log  {  Serious about programming in Tcl ?}                             {hype heading}
    log  {  }                                                                hype
    log* {  Get }                                                            hype
    log*        {ActiveTcl Pro Studio}                                      {hype emphasis}
    log                             {, the complete Tcl development bundle.} hype
    log* {  }                                                                hype
    log*   {ActiveTcl Pro Studio}                                           {hype emphasis}
    log                        { combines everything you need in a single, } hype
    log  {  high-value bundle, including advanced Tcl programming}           hype
    log  {  technologies and in-depth resource information.}                 hype
    log  {  }                                                                hype
    log  {  Find out more:}                                                  hype
    log* {    }                                                              hype
    log*      { http://www.ActiveState.com/Tcl }                            {hype url}
    log                                               { }                    hype
    log  {  }                                                                hype
}

proc teapot_message {} {
    if {![string equal $::AT(MODE) normal]} return
    if {[package vcompare $::AT(tcl_VERSION) 8.5] < 0} return

    log  { }
    log  { }                                                          {note}
    log* "  Note that ActiveTcl $::AT(VERSION) is a "                 {note}
    log*                                            {trimmed down}    {note emphasis}
    log                                                          { distribution} {note}
    #log                { providing only the most important packages.} {note}
    log  {  providing only the most important packages. All packages}    {note}
    log  {  not found in the distribution can be installed by using}     {note}
    log* {  the }                                                        {note}
    log*       {teacup client}                                           {note emphasis}
    log*                    { to the }                                   {note}
    log*                            {TEApot Package Management}          {note emphasis}
    log                                                      { however.} {note}
    log  { }                                                             {note}
    log  {  Further note that the documentation was not trimmed, and}    {note}
    log* {  contains the documentation of }                              {note}
    log*                                 {all}                           {note emphasis}
    log                                     { packages, even those not}  {note}
    log  {  installed by the distribution.}                              {note}
    log  { }                                                             {note}
    log  { }
}

# df -- Determine available free space for a directory.
#
# Returns value in KB as a double.
#
proc df {{dir .}} {
    # Force the result to be a double, to allow 51 bits of int precision
    # and not error on large disk allocations
    switch $::tcl_platform(os) {
	Darwin - FreeBSD - Linux - OSF1 - SunOS {
	    # Bugfix: If the path in an entry is too long SunOS splits
	    # the entry into two lines, one containing the path, the
	    # other containing the remainder of the information. This
	    # means that every in the last line moves one index to the
	    # front and index 3 refers to the percentage instead of
	    # the available space. Our fix is to count the elements
	    # from the end, as the available space is the third one
	    # from the end, and this index does not change when an
	    # entry is split into two lines.

	    return [lindex [lindex [split [exec df -k $dir] \n] end] end-2].0
	}
	HP-UX {
	    # Bugfix. See Linux etc. above. The same can happen on HP-UX.
	    # Thanks to Nick Fines <nick.fines@marconi.com>

	    return [lindex [lindex [split [exec bdf   $dir] \n] end] end-2].0
	}
	AIX {
	    return [lindex [lindex [split [exec df -k $dir] \n] end] 2].0
	}
	{Windows NT} - {Windows 95} {
	    catch {
		set dir [file nativename $dir]
		set res [eval exec [auto_execok dir] [list $dir]]
		set line [lindex [split $res "\n"] end]
		if {[regexp -nocase {([0-9,\.]+)\s+(bytes|KB|MB)\s+} \
			$line -> size type]} {
		    set size [string map {, {}} $size]
		    switch $type {
			MB - mb { return [expr {$size * 1000.0}] }
			KB - kb { return [expr {double($size)}] }
			BYTES - bytes {
			    if {[string match {*.*} $size]} {
				return [expr {$size / 1024.0}]
			    } else {
				return [expr {double($size) / 1024.0}]
			    }
			}
		    }
		} else {
		    error "Unable to extract free space for $dir from output of 'dir'"
		}
	    }
	    # Some error occured, assume we have at least 100MB
	    return 100000.0
	}
	default {error "Unable to get disk free space on $::tcl_platform(os)"}
   }
}

# Determine the used space for a directory.
#
proc du {{dir .}} {
    switch -glob $::tcl_platform(os) {
	Darwin - FreeBSD - Linux - OSF1 - SunOS - HP-UX - AIX {
	    return [lindex [exec du -sk $dir] 0]
	}
	"Windows*" {
	    # Just say we need 25MB
	    return 25000
	}
	default {error "Unable get disk usage on $::tcl_platform(os)"}
    }
}

# Check if enough space for $src dir exists in $dest dir
# Returns 0 on OK, otherwise space required
#
proc need_space {src dest} {
    set need [du $src]
    set have [df $dest]
    return [expr {($need > $have) ? $need : 0}]
}

proc check_installdir_basic {installdir} {
    global INSTALL_MODE tcl_platform

    if {[string equal {} $installdir]} {
	set installdir [default_installdir]
	set errmsg "Please choose a non-empty directory name."
	return [list $installdir $errmsg]
    }

    # Change the directory into an absolute path for better comparison
    # later.

    set installdir [file join [pwd] $installdir]

    # Check installation dir for existence and correct type.

    if {[file exists $installdir] && ![file isdirectory $installdir]} {
	set errmsg "The chosen path \"$installdir\" is not a directory.\
		\nPlease choose a directory."
	return [list $installdir $errmsg]
    } elseif {[string equal $installdir $::SCRIPT_DIR]} {
	set errmsg "You are trying to install $::AT(dPGROUP) over the directory\
		\ncontaining the distributed files. This is not allowed.\
		\nPlease choose a different directory."
	return [list $installdir $errmsg]
    } elseif {[inside_distribution $installdir]} {
	set errmsg "You are trying to install $::AT(dPGROUP) inside of the directory\
		\ncontaining the distributed files. This is not allowed.\
		\nPlease choose a different directory."
	return [list $installdir $errmsg]
    }

    return [list $installdir {}]
}

proc check_installdir_post {installdir} {
    set errmsg ""
    if {[catch {
	if {![file exists $installdir]} { file mkdir $installdir }
	set pwd [pwd]
	cd $installdir
	cd $pwd
    } err]} {
	set errmsg "Invalid directory choice: $err.\
		\nPlease choose another directory."
    } elseif {![file writable $installdir]} {
	set errmsg "Cannot write to directory \"$installdir\".\
		\nPlease choose another directory."
    } elseif {[set need [need_space $::SCRIPT_DIR $installdir]]} {
	# Check of the available space for the chosen directory against
	# our requirements. Instead of using a fixed requirement we
	# query the unpacked archive for its size.

	# The check is conservative as it counts the demos as part of
	# the requirements even if they will be installed in a
	# different location.
	set errmsg "The chosen path \"$installdir\" does not provide\
		\nenough space for the installation of $::AT(dPGROUP).\
		\n\nPlease choose a directory with at least ${need}K\
		available space.\n"
    }

    return [list $installdir $errmsg]
}


proc check_installdir {installdir} {
    set ::AT(InstVersion) ""
    set errmsg ""

    foreach {dir msg} [check_installdir_basic $installdir] break
    if {$msg != {}} {return [list $dir $msg]}

    set installdir $dir

    # Check installation dir for existence and correct type.

    if {[file exists $installdir]} {
	# Code executed for ActiveTcl and full Tcl/DevKit installation.

	if {[string equal $::AT(MODE) normal]} {
	    # Look for an existing installation in the chosen directory.
	    # Issue a warning even if overwriting is allowed.

	    # NOTE: We allow the installation of ActiveTcl even if the
	    # chosen directory contains a Tcl/DevKit installation. This
	    # allows the user/customer to replace the ActiveTcl part
	    # of a Tcl/DevKit installation with an updated version.

	    if {
		[contains_activetcl_installation $installdir] &&
		![overinstall_allowed]
	    } {
		set errmsg "The chosen directory contains a different major installation\
			\nof $::AT(dPGROUP) ($::AT(InstVersion)).\
			\nPlease choose a different directory."
		return [list $installdir $errmsg]
	    }
	} else {
	    # Look for an existing installation in the chosen directory.
	    # Issue a warning even if overwriting is allowed.

	    # NOTE: We allow the installation of a full Tcl/DevKit even if the
	    # chosen directory contains an ActiveTcl installation. This
	    # allows the user/customer to upgrade to Tcl/DevKit while
	    # simultaneously (possibly) upgrading the ActiveTcl
	    # part. There is no need to force the user to use the
	    # Upgrade Kit in this case. It is less efficient regarding
	    # download time and required disk space, but well, that is
	    # the users choice.

	    if {
		[contains_tdk_installation $installdir] &&
		![overinstall_allowed]
	    } {
		set errmsg "The chosen directory contains a different major installation\
			\nof $::AT(dPGROUP) ($::AT(InstVersion)).\
			\nPlease choose a different directory."
		return [list $installdir $errmsg]
	    }
	}
    }

    return [check_installdir_post $installdir]
}

proc inside_distribution {installdir} {
    set plen    [llength [file split $::SCRIPT_DIR]]
    set iprefix [eval [list file join] [lrange [file split $installdir] \
	    0 [incr plen -1]]]

    return [string equal $iprefix $::SCRIPT_DIR]
}

proc check_demodir {demodir} {
    set errmsg  ""
    if {$demodir == ""} {
	set demodir [default_demodir $::INSTALL_DIR]
	set errmsg "Please choose a non-empty directory name."
    } elseif {[file exists $demodir] && ![file isdirectory $demodir]} {
	# Check of installation dir for existence and correct type.
	set errmsg "The chosen path \"$demodir\" is not a directory.\
		\nPlease choose a directory."
    } elseif {[catch {
	set pwd [pwd]
	if {![file exists $demodir]} { file mkdir $demodir }
	cd $demodir
	cd $pwd
    } err]} {
	set errmsg "Invalid directory choice: $err.\
		\nPlease choose another directory."
    } elseif {![file writable $demodir]} {
	set errmsg "Cannot write to directory \"$demodir\".\
		\nPlease choose another directory."
    } elseif {[set need [need_space [file join $::SCRIPT_DIR pdemos] $demodir]]} {
	# Check the available space for the chosen directory against
	# our requirements. Instead of using a fixed requirement we
	# query the unpacked archive for its size.
	set errmsg "The chosen path \"$demodir\" does not provide\
		\nenough space for the installation of demo scripts.\
		\n\nPlease choose a directory with at least ${need}K\
		available space.\n"
    }

    return [list $demodir $errmsg]
}

proc check_runtimedir {runtimedir} {
    set errmsg  ""
    if {$runtimedir == ""} {
	set runtimedir $::INSTALL_DIR
	set errmsg "Please choose a non-empty directory name."
    } elseif {[file exists $runtimedir] && ![file isdirectory $runtimedir]} {
	# Check of installation dir for correct type, if existing
	set errmsg "The chosen path \"$runtimedir\" is not a directory.\
		\nPlease choose a directory."
    }
    # Note: The runtime dir may not exist. This is ok, as it is the
    # directory the final user will see from his perspective of the
    # filesystem. The installer may not see it at all.

    return [list $runtimedir $errmsg]
}

proc default_installdir {} {
    if {
	[string equal $::AT(MODE) normal]
    } {
	# ActiveTcl

	if {[string equal "windows" $::tcl_platform(platform)]} {

	    if {[info exists env(SYSTEMDRIVE)]} {
		return [file nativename "$env(SYSTEMDRIVE)\\Tcl"]
	    } else {
		return [file nativename "C:\\Tcl"]
	    }
	} else {
	    set major [join [lrange [split $::AT(VERSION) .] 0 1] .]
	    return "/opt/ActiveTcl-$major"
	}
    } elseif {[string equal $::AT(MODE) pro]} {
	# TclDevKit

	if {[string equal "windows" $::tcl_platform(platform)]} {
	    return [file nativename "C:\\TclDevKit"]
	} else {
	    set major [join [lrange [split $::AT(VERSION) .] 0 1] .]
	    return "/opt/TclDevKit-$major"
	}
    } else {
	return -code error "PANIC. Invalid distribution: $::AT(MODE)"
    }
}

proc default_demodir {installdir} {
    return [file nativename [file join $installdir demos]]
}

proc install_log {msg} {
    if {[info exists ::INSTALL_DIR]} {
	if {![info exists ::INSTALL_LOG]} {
	    global INSTALL_MODE

	    # The first call to log an action taken by the installer
	    # creates and initializes the directory to contain the log
	    # and the uninstaller application. Reduced set of images,
	    # as the uninstaller is wrapped, i.e. it carries them
	    # within itself. Only the image needed for the windows
	    # shortcuts is put into the relevant location.

	    set logdir [file dirname $::INSTALL_DIR/$::AT(LOG)]

	    file mkdir $logdir $logdir/install_images
	    file copy -force \
		    $::SCRIPT_DIR/$::AT(UNINSTALLER) \
		    $logdir/$::AT(UNINSTALLER)

	    file copy -force \
		    $::SCRIPT_DIR/install_images/uninstall.ico \
		    $logdir/install_images/uninstall.ico

	    set ::INSTALL_LOG [open $::INSTALL_DIR/$::AT(LOG) w]

	    puts $::INSTALL_LOG "# -*- tcl -*- $::AT(LOG) v1.0\
		    \n# DO NOT DELETE THIS FILE.\
		    \n# It is needed by the uninstaller to clean up properly.\
		    \n# Installed on [clock format [clock seconds]]\n"

	    # Remember information about the type of installation
	    # (common / personal, admin rights).

	    puts $::INSTALL_LOG [list INSTALLTYPE $INSTALL_MODE]

	    # Fake the create entries for the directory containing the log file.
	    set d $logdir
	    set ndir [file nativename $::INSTALL_DIR]
	    while {![string equal [file nativename $d] $ndir]} {
		puts $::INSTALL_LOG [list MKDIR $d]
		set d [file dirname $d]
	    }
	    puts $::INSTALL_LOG [list FILECOPY $::INSTALL_DIR/$::AT(LOG)]
	    puts $::INSTALL_LOG [list FILECOPY $logdir/$::AT(UNINSTALLER)]

	    puts $::INSTALL_LOG [list MKDIR    $logdir/install_images]
	    puts $::INSTALL_LOG [list FILECOPY $logdir/install_images/uninstall.ico]
	}
	puts $::INSTALL_LOG $msg
    }
}

proc install_all {srcdir installdir {srcrelative {}}} {
    # Assumes existence of a [log] command.

    if {$srcrelative eq "lib/teapot"} {
	switch -exact -- $::INSTALL_SEED {
	    existing {
		# This branch is run not only for seed and repo in
		# conflict, but also if there is no conflict at all.
		if {[file exists $installdir]} {
		    log "\tA repository already exists, and you have chosen to keep it."  url
		    log "\tSkipping the seed repository in directory" url
		    log "\t\t[file nativename $srcdir]" url
		    return
		}
	    }
	    seed {
		log "\tA repository already exists, and you have chosen to drop it"  url
		log "\tin favor of the seed repository. Renaming it to teapot-old."  url
		log "\tAn existing teapot-old directory will be removed."            url

		while {[catch {
		    file delete -force ${installdir}-old
		    file rename $installdir ${installdir}-old
		} msg]} {
		    switch -exact [abortretryignore "Saving the repository failed.\nThe operating system reported:\n\t$msg\n"] {
			abort  { exit 1 }
			retry  continue
			ignore break
		    }
		}
	    }
	    merge {
		log "\tA repository already exists, and you have chosen to merge"  url
		log "\tthe seed repository into it. This will happen later as part" url
		log "\tof the general repository setup." url
		return
	    }
	}
    }

    if {![file exists $installdir]} {
	log "\tCreating directory [file nativename $installdir] ..."
	file mkdir $installdir
    }

    # Note in the install log that we "made" the directory in either case,
    # so it will try to delete that directory if it is empty.

    install_log [list MKDIR $installdir]

    foreach f [glob -nocomplain -directory $srcdir -types {f l} *] {
	#log "\t[file join $installdir [file tail $f]] ..."
	if {[catch {
	    file copy -force $f $installdir
	    install_log [list FILECOPY \
		    [file join $installdir [file tail $f]]]
	} err]} { log $err error }
    }

    foreach dir [glob -nocomplain -directory $srcdir -type d *] {
	set tail [file tail $dir]
	set inst [file join $installdir $tail]

	install_all $dir $inst [file join $srcrelative $tail]
    }
}

proc patch_txr {f atdemodir installdir} {
    set map [list \
		 @atdemodir@  [file nativename $atdemodir] \
		 @installdir@ [file nativename $installdir]]

    mk::file open db $f
    # Keep in sync with devkit/app/xref/2mk.tcl
    mk::view layout db.file { md5:S where_str:S rowid:I }
    # Loop over the files and map the placeholders
    mk::loop c db.file {
	mk::set $c where_str [file nativename [string map $map [mk::get $c where_str]]]
    }
    mk::file commit db
    mk::file close  db
    return
}

proc patch_tools {installdir runtimedir} {
    # We can savely skip this action if the distribution does not
    # contain starkits, or doesn't need patching (platform).
    # - Rewrite - Normal and without Starkits have 'Tkcon' to be
    #   patched equivalently

    if {[string equal "windows" $::tcl_platform(platform)]} {return}

    if {![string equal $::AT(MODE) pro]} {
	# Normal/Regular distribution. Unix. Patch Tkcon with the path
	# to the installed wish. Do not do this for TDK. TDK does not
	# come with Tkcon, and would touch an already installed one if
	# TDK is put into the same directory as AT. Mangling the
	# installed tkcon, causing it to not work anymore.

	# The same reasoning applies to AT 8.5. It has no Tkcon either
	# and would mangle an installed one as well.
	# This part has to be checked, does it work with overinstall ?
	#
	#if {[package vsatisfies $::AT(tcl_VERSION) 8.5]} return

	log "-   Patching AT applications" emphasis

	set wfiles [list]
	foreach {app word} {
	    tkcon      wish
	    tclvfse    tclsh
	    guibuilder tclsh
	} {
	    set f [file join $installdir bin $app]
	    if {![file exists $f]} continue
	    log "*  [file nativename $f] ... (wish)"
	    patch_kit $f $word [file join $runtimedir bin wish$::AT(tcl_VERSION)]
	}

	if {[package vsatisfies $::AT(tcl_VERSION) 8.5]} return
	# 8.4 only

	foreach f {
	    dtplite page tcldocstrip
	} {
	    set fx [file join $installdir bin $f]
	    if {![file exists $fx]} continue
	    log "*  [file nativename $fx] ... (tclsh)"
	    patch_kit $fx tclsh [file join $runtimedir bin tclsh]
	}
    } else {
	# Pro distribution aka TDK. Unix. Patch the custom TDK basekit
	# into the tool starkits for proper execution.

	set files {}
	# tclvfse
	foreach t {
	    tclapp tclchecker tclcompiler tcldebugger
	    tclinspector tclpe tclxref
	} {
	    set t [file join $installdir bin $t]
	    if {[file exists $t]} {
		lappend files $t
	    }
	}

	log "-   Patching TDK Tools" emphasis

	foreach f $files {
	    log "*  [file nativename $f] ... (tdkbase)"
	    patch_kit $f tclsh [file join $runtimedir bin tdkbase]
	}
    }
    return
}

proc patch_kit {kit word new} {
    # Replaces the first occurence of 'word' with 'new', in the header
    # (!) of the specified star'kit'.
    #
    # Unix-only!

    set tmpfile $kit.[pid]
    set out [open $tmpfile w]
    set in  [open $kit r]

    while {1} {
	gets $in line
	if {[set start [string first ${word} $line]] >= 0} {
	    set end [expr {$start + [string length $word]}]
	    incr start -1
	    puts -nonewline $out [string range $line 0 $start]

	    # Bugzilla 23548. Quote to handle spaces in the path.
	    puts -nonewline $out \"$new\"

	    puts            $out [string range $line $end end]
	    break
	}
	puts $out $line
    }
    fconfigure $out -encoding binary -translation binary
    fconfigure $in  -encoding binary -translation binary
    fcopy $in $out
    close $in
    close $out
    file copy -force $tmpfile $kit
    file delete $tmpfile

    if {[string equal $::tcl_platform(platform) windows]} {
	file attributes $kit -readonly 0
    } else {
	file attributes $kit -permissions 0755
    }
    return
}



proc patch_tut {} {
    if {[string equal "windows" $::tcl_platform(platform)]} {
	set exeext .exe
	set dllext dll
	set tclext .tcl
    } else {
	set tclext ""
	set exeext ""
	set dllext so
	if {[string equal $::tcl_platform(os) HP-UX]} {
	    set dllext sl
	}
    }
    if {$::AT(DEBUG)} {
	set exeext g$exeext
    }

    # Bug 69643. For TDK the threaded status of itself is _irrelevant_,
    # we have to look at the status of AT to use the proper basekits.
    # Determine this indirectly, through querying the filesystem
    # (glob).

    if {[string equal $::AT(MODE) pro]} {

	set base     [file join $::ATinstalldir bin]
	set basekits [glob -nocomplain -directory $base base-tcl*]

	if {[lsearch -glob $basekits *-thread-*] < 0} {
	    set threadext {}
	} else {
	    set threadext {thread-}
	}
    } elseif {$::AT(THREADED)} {
	set threadext {thread-}
    } else {
	set threadext {}
    }

    if {
	![string equal $::AT(MODE) normal]
    } {
	## Patch project files in debugger tutorial code.
	## Pro and lite only.

	## In lite the demos are not present, hence nocomplain to prevent
	## the installer from balking.

	set demodir [file join $::DEMO_DIR TclDevKit]

	set files [glob -nocomplain -directory $demodir */*/*.tpj]

	log "-   Patching TDK project files" emphasis

	foreach f $files {
	    log "*  [file nativename $f] ..."
	    patch_file_simple $f @demodir@    $::DEMO_DIR
	    patch_file_simple $f @installdir@ $::RUNTIME_DIR
	    patch_file_simple $f @exeext@     $exeext
	    patch_file_simple $f @arch@       $::AT(ARCH)
	    patch_file_simple $f @thread@     $threadext
	    patch_file_simple $f @tclext@     $tclext
	    patch_file_simple $f @dllext@     $dllext

	    # Some parts of some demos rely on the existence of an
	    # ActiveTcl installation, using files it provides
	    # (basekits, AT demos). See [xx] too for definition.

	    patch_file_simple $f @atdemodir@    $::ATdemodir
	    patch_file_simple $f @atinstalldir@ $::ATinstalldir
	}

	# Patch example cross-reference databases.

	set files [glob -nocomplain -directory $demodir */*.txr]

	log "-   Patching TDK cross-reference databases" emphasis

	foreach f $files {
	    log "*  [file nativename $f] ..."

	    patch_txr $f $::ATdemodir $::RUNTIME_DIR
	}
    }

    if 0 {
	if {
	    [string equal $::AT(MODE) pro] &&
	    [string equal $::AT(tcl_VERSION) 8.4]
	} {
	    log "-   Patching Expect demos" emphasis

	    set demodir $::DEMO_DIR
	    set files [glob -nocomplain -directory [file join $demodir TclDevKit TclApp Expect] *.tpj]

	    foreach f $files {
		log "*  [file nativename $f] ..."
		patch_file_simple $f @demodir@      $::DEMO_DIR
		patch_file_simple $f @installdir@   $::RUNTIME_DIR
		patch_file_simple $f @atdemodir@    $::ATdemodir
		patch_file_simple $f @atinstalldir@ $::ATinstalldir
		patch_file_simple $f @exeext@       $exeext
		patch_file_simple $f @arch@         $::AT(ARCH)
		patch_file_simple $f @thread@       $threadext
		patch_file_simple $f @tclext@       $tclext
		patch_file_simple $f @dllext@       $dllext
	    }
	}
    }
    return
}


proc placeholder {} {
    # see also make-active.tcl
    return @[string repeat _ 400]@
}

proc patch_shells {srcdir installdir runtimedir} {
    # No patching of anything for the lite applications.
    # WRONG. This was introduced in rev 40, but is wrong. We have to
    # patch the distributed proshells and wrapstubs like for the pro
    # distribution. Conditialize each section separately.
    ## if {[string equal $::AT(MODE) lite]} {return} ##

    set longfakepath [placeholder]

    # Check string sizes before trying to patch.

    if {[string length $installdir] > [string length $longfakepath]} {
	log "Can't patch shells and libraries, new path to long"
	return
    }

    set libdir    [file join $runtimedir lib]
    set libdirtcl [file join $runtimedir lib tcl$::AT(tcl_VERSION)]

    if {![string equal $::AT(MODE) pro]} {
	if {$::AT(DEBUG)} { set dbg "g" } else { set dbg "" }

	set files [list]

	# Always patch the versioned shells
	set patterns [list tclsh${dbg}$::AT(tcl_VERSION) \
			   wish${dbg}$::AT(tcl_VERSION)  \
			   tclsh$::AT(tcl_VERSION)${dbg} \
			   wish$::AT(tcl_VERSION)${dbg}]

	# Disabled. Have unversioned always now. Patch unversioned shells only for 8.4
	if {1||![package vsatisfies $::AT(tcl_VERSION) 8.5]} {
	    lappend patterns tclsh${dbg} wish${dbg}
	}

	foreach f $patterns {
	    set ff [file join $installdir bin $f]
	    if {[file exists $ff]} {
		lappend files $ff
	    }
	}

	foreach f $files {
	    log "* [file nativename $f] ..."
	    patch_file $f $longfakepath $libdir
	}
    }

    set files [list]
    if {![string equal $::AT(MODE) pro]} {
	# Patch the Tcl shared and static libraries
	# Note: Only for the version we are currently installing.
	set files [glob -nocomplain -directory [file join $installdir lib] \
		       "libtcl$::AT(tcl_VERSION)*.{so,sl,a}"]
    }

    ### The patching of the main tcl library is different between 8.3 and 8.4,
    ### and 8.5.  Also special for HP, AIX, and Solaris/IA32 (x86)

    if {[package vsatisfies $::AT(tcl_VERSION) 8.5]} {
	# Patching 8.5

	if {"$::tcl_platform(os)" == "AIX"} {
	    # AIX is special. The fourth item needs the special AIX handling
	    # for path (keep :). The first has to be regular for proper
	    # auto_path. The other two can be either/or with affecting things.

	    foreach f $files {
		log "*  AIX 8.5 [file nativename $f] ..."
		patch_file $f $longfakepath \
			$libdir \
			$libdirtcl \
			$libdirtcl \
			$libdir
	    }
	} else {
	    foreach f $files {
		log "*  8.5 [file nativename $f] ..."
		patch_file $f $longfakepath \
			$libdir \
			$libdirtcl \
			$libdirtcl \
			$libdir
	    }
	}
    } elseif {("$::tcl_platform(os)" == "HP-UX")} {
	foreach f $files {
	    log "*  HP [file nativename $f] ..."
	    patch_file $f $longfakepath \
		    $libdir \
		    $libdirtcl \
		    $libdirtcl
	}
    } elseif {
	("$::AT(ARCH)" == "solaris-ix86") &&
	[package vsatisfies $::AT(tcl_VERSION) 8.4]
    } {
	# Patching 8.4 solaris intel
	foreach f $files {
	    log "*  SOL 8.4 x86 [file nativename $f] ..."
	    patch_file $f $longfakepath \
		$libdir \
		$libdir \
		$libdirtcl
	}
    } elseif {("$::tcl_platform(os)" == "AIX") ||
	      [package vsatisfies $::AT(tcl_VERSION) 8.4]} {
	# Patching 8.4 (or AIX)
	foreach f $files {
	    log "*  AIX|8.4 [file nativename $f] ..."
	    patch_file $f $longfakepath \
		$libdir \
		$libdirtcl \
		$libdirtcl
	}
    } else {
	# Patching 8.3
	foreach f $files {
	    log "*  [file nativename $f] ..."
	    patch_file $f $longfakepath \
		    $libdirtcl \
		    $libdir
	}
    }

    if {![string equal $::AT(MODE) pro]} {
	# Patch the Tk shared and static libraries
	set files [glob -nocomplain -directory [file join $installdir lib] \
		       "libtk$::AT(tcl_VERSION)*.{so,sl,a}"]

	foreach f $files {
	    log "*  [file nativename $f] ..."
	    patch_file $f $longfakepath $libdir
	}
    }

    ### ######### ###########################
    ## Expect: Now patch the expect libraries as well, they also
    ## contain a path to their own location, exposed at the tcl level
    ## through the variable --> 'expect_library'.
    ## This is 8.4 only, as 8.5 does not come with Expect!.
    
    set el [glob -nocomplain -directory [file join $installdir lib] expect\[0-9\]*]
    if {[llength $el]} {
	# Sort in decreasing order to get the newest expect for
	# patching in case when multiple expect dirs are
	# found. This can happen if a new ActiveTcl distribution
	# is installed over an excisting older one containing an
	# older expect.
	set el [lsort -dictionary -decreasing $el]

	regsub -- {^expect} [file tail [set el [lindex $el 0]]] {} eversion

	set files [glob -nocomplain -directory $el {libexpect[0-9]*.{so,a,sl}}]

	# When patching the library the first placeholder is for
	# -rpath, the lib dir, and the second one refers to the
	# directory containing the expect tcl library.

	if {$::tcl_platform(os) == "AIX"} {
	    # On AIX we see 4 placeholders, and we must not modify
	    # the second to last, only the first.
	    foreach f $files {
		log "*  [file nativename $f] ..."
		patch_file_first $f $longfakepath \
		    [file join $runtimedir lib expect$eversion]
	    }
	} else {
	    foreach f $files {
		log "*  [file nativename $f] ..."
		patch_file $f $longfakepath \
		    [file join $runtimedir lib expect$eversion]
	    }
	}

	# We patch the expect shells as well, if they are present.
	# Usually they are not, but customized distributions may
	# contain them.

	set files [glob -nocomplain -directory [file join $installdir bin]\
		       expect*]

	# When patching the shells we have two placeholders, which
	# are for -rpath, and have to contains the lib dir, and
	# package dir. The placeholders are in one string, and we
	# do not wish to cut after the first, so we replace both
	# in one. Note!! The second place holder already ends in
	# /expectX.Y. Using expectX.y in the first placeholder is
	# therefore bad! Leave it out. The first is lib dir, the
	# second + trailer is package dir.

	foreach f $files {
	    log "*  [file nativename $f] ..."
	    patch_file $f \
		${longfakepath}:${longfakepath} \
		[file join $runtimedir lib]:$libdir
	}
    }

    ## Expect: Patching done ...
    ### ######### ###########################

    if {1} {
	### ######### ###########################
	## Scotty: Patch the Tnm library to contain the paths to its helper demon
	## ......: applications. Otherwise a path is required.

	set el [glob -nocomplain -directory [file join $installdir lib] \
		    {tnm[0-9]*}]
	if {[llength $el]} {
	    set tnm [lindex $el 0]

	    set files [glob -nocomplain -directory $tnm \
			   {libtnm[0-9]*.{so,sl,a}}]
	    foreach f $files {
		log "*  [file nativename $f] ..."
		patch_file $f icmp$longfakepath [file join $runtimedir bin/nmicmpd]
		patch_file $f trap$longfakepath [file join $runtimedir bin/nmtrapd]
	    }
	}
    }

    ## Scotty/Tnm: Patching done ...
    ### ######### ###########################

    ## Patch the *Config.sh files too
    if {![string equal $::AT(MODE) pro]} {

	log "Patching configuration scripts ($::AT(tcl_VERSION))"

	set files [glob -nocomplain -directory [file join $installdir lib] \
		       {t{cl,k}Config*.sh} \
		       "t{cl,k}$::AT(tcl_VERSION)/t{cl,k}Config*.sh"]
	foreach f $files {
	    log "*  [file nativename $f] ..."
	    patch_text_file $f $longfakepath $runtimedir
	}
    }
    return
}

proc patch_text_file {file key value} {
    ## patch a text file mapping $key to $value
    ## log "   - [list $key | $value]"

    if {[string equal $::tcl_platform(platform) windows]} {
	file attributes $file -readonly 0
    } else {
	file attributes $file -permissions u+w
    }

    set in  [open $file r]
    set out [open $file.tmp w]

    ## If the key contains backslashes (windows paths) match them
    ## literally!

    set key [string map [list \\ \\\\] $key]

    while {[gets $in line] >= 0} {
	regsub -all -- $key $line $value line
	puts $out $line
    }
    close $in
    close $out
    file rename -force $file.tmp $file
}

proc patch_file_simple {file key value} {
    set data [read [set fh [open $file r]]]
    close $fh

    set data [string map [list $key $value] $data]

    set              fh [open $file w]
    puts -nonewline $fh $data
    close           $fh
    return
}

proc patch_pad {data key value {all 0}} {
    # Always search for key until next \0 (null), replace key with
    # correct value and pad with \0 to previous length
    if {$all} { set opt -all } else { set opt -indices }
    foreach indx [regexp -inline -indices $opt -- "$key.*?\0" $data] {
	foreach {start end} $indx { break }
	incr end -1 ; # keep null
	set match [string range $data $start $end]
	set repl [string map [list $key $value] $match]
	set mlen [string length $match]
	set rlen [string length $repl]
	set pad "\0"
	if {$::tcl_platform(os) == "AIX" && [string match *:* $match]} {
	    # if we are on AIX and have a multi-part path, then this needs
	    # special case handling to pad with : instead of null
	    set pad ":"
	}
	append repl [string repeat $pad [expr {$mlen-$rlen}]]
	set data [string replace $data $start $end $repl]
    }
    return $data
}

proc patch_file {file key args} {
    if {![llength $args]} {
	return -code error "Need at least one patch value"
    }
    foreach value $args {
	if {[string length $value] > [string length $key]} {
	    log "Can't patch \"$file\", value too long"
	    return
	}
    }

    # Patch inplace ... (and make sure that we can!)
    # Read everything into memory, patch it, seek to beginning of file
    # and write the modifications back.

    if {[string equal $::tcl_platform(platform) windows]} {
	file attributes $file -readonly 0
    } else {
	file attributes $file -permissions u+w
    }

    set fh [open $file r+]
    fconfigure $fh -encoding binary -translation binary

    set data [read $fh]

    # More than one value. Patch one by one ...
    set lastval [lindex $args end]
    foreach value [lrange $args 0 end-1] {
	set data [patch_pad $data $key $value]
    }
    set data [patch_pad $data $key $lastval 1] ; # do all the rest

    seek $fh 0 start
    puts -nonewline $fh $data
    close $fh

    if {[string equal $::tcl_platform(platform) windows]} {
	file attributes $file -readonly 1
    } else {
	file attributes $file -permissions og-w
    }
    return
}

proc patch_file_first {file key value} {
    if {[string length $value] > [string length $key]} {
	log "Can't patch \"$file\", value too long"
	return
    }

    # Patch inplace ... (and make sure that we can!)
    # Read everything into memory, patch it, seek to beginning of file
    # and write the modifications back.

    if {[string equal $::tcl_platform(platform) windows]} {
	file attributes $file -readonly 0
    } else {
	file attributes $file -permissions u+w
    }

    set fh [open $file r+]
    fconfigure $fh -encoding binary -translation binary

    set data [patch_pad [read $fh] $key $value]

    seek $fh 0 start
    puts -nonewline $fh $data
    close $fh

    if {[string equal $::tcl_platform(platform) windows]} {
	file attributes $file -readonly 1
    } else {
	file attributes $file -permissions og-w
    }
    return
}

proc set_key {key value data} {
    # make sure to log when we set a key
    r_set $key $value $data
    install_log [list REGKEY $key $value $data]
}

proc r_set {key value data {mod {}}} {
    if {$mod == {}} {
	log_debug "\treg set \[$key, $value\] = ($data)"
	registry set $key $value $data
    } else {
	log_debug "\treg set ($mod) \[$key, $value\] = ($data)"
	registry set $key $value $data $mod
    }
}

proc r_get {key value} {
    log_debug "\treg get \[$key, $value\] = "

    set res [registry get $key $value]

    log_debug "\t\t$res"
    return $res
}

proc get_editor_app {} {
    set write [file nativename [auto_execok write]]
    if {$write == ""} {
	set write notepad.exe
    }
    return $write
}

proc get_edit_assoc {assoc ext} {
    # If an existing Edit association is there, don't change it,
    # as it may be a special one set up by the user.

    if {[catch {
	r_get HKEY_CLASSES_ROOT\\$assoc\\shell\\edit\\command {}
    } write]} {
	set write "\"[get_editor_app]\" \"%1\""
    } else {
	log "\tPreserving $ext Edit association \"$write\""
    }

    return $write
}


proc set_tcl_association {installdir ndir bindir} {
    ## Create the file association
    ## Not required for lite, we expect an existing ActiveTcl install
    ## and thus assume the existence of the association

    if {[string equal $::AT(MODE) pro]}        {return}

    set assoc ActiveTclScript

    if {[catch {
	## Make sure the user requested that we install this association.
	if {![info exists ::INSTALL_REG_TCL] || $::INSTALL_REG_TCL} {
	    log "\tAdding recognition of .tcl and .tbc scripts ..."

	    # If an existing Edit association is there, don't change it,
	    # as it may be a special one set up by the user.
	    set write [get_edit_assoc $assoc .tcl]
	    set assocList {}
	    lappend assocList \
		$assoc			{} "ActiveTcl Script" \
		$assoc\\DefaultIcon	{} "$bindir\\tk$::AT(ND_VERG).dll,0" \
		$assoc\\shell\\edit	{} &Edit \
		$assoc\\shell\\edit\\command {} $write \
		$assoc\\shell\\open\\command {} "\"$bindir\\wish$::AT(ND_VERG).exe\" \"%1\" %*"
	    lappend assocList .tcl {} $assoc
	    # Only allow .tbc if .tcl is also associated
	    if {![info exists ::INSTALL_REG_TBC] || $::INSTALL_REG_TBC} {
		lappend assocList .tbc {} $assoc
	    }
	    foreach {key valueName data} $assocList {
		set_key [class_root]\\$key $valueName $data
	    }
	}
    } err]} {
	log "ERROR: $err" error
    }
    return
}


proc set_association {
    assoc title
    application arguments
    extensions
    iconsource
} {
    ## Create file associations
    ## Assume that the application is a list. One element => True application
    ## More elements => We did resolve chaining of open commands
    #
    ## Required for starkit based launcher. windows does not search
    ## for an open command for the open application of a file.

    set bindir [file nativename $::INSTALL_DIR]\\bin

    if {[catch {
	#log "___ $application"

	if {[llength $application] > 1} {
	    # Chained app, quote all parts.
	    set application "\"$bindir\\[join $application "\" \"$bindir\\" ]\""
	} elseif {[llength $application] == 1} {
	    # Normal app
	    set application "\"$bindir\\[lindex $application 0]\""
	} else {
	    # error, illegal call
	    error "internal error in call to 'set_association'"
	}

	# Now tack on all the arguments which do not require knowledge
	# of the bin directory, just regular quoting.
	foreach a $arguments {append application " \"$a\""}

	#log "___ $application"

	set write [get_edit_assoc $assoc $assoc]

	foreach {key valueName data} [list \
		$assoc			     {} "$title" \
		$assoc\\DefaultIcon	     {} $iconsource \
		$assoc\\shell\\edit	     {} &Edit \
		$assoc\\shell\\edit\\command {} $write \
		$assoc\\shell\\open\\command {} "$application \"%1\" %*" \
		] {
	    set_key [class_root]\\$key $valueName $data
	}

	foreach ext $extensions {
	    log "\tAdding recognition of $ext files ..."

	    foreach {key valueName data} [list \
		    $ext {} "$assoc" \
		    ] {
		set_key [class_root]\\$key $valueName $data
	    }
	}
    } err]} {
	log "ERROR: $err" error
    }
    return
}




proc set_product {installdir ndir bindir} {
    ## Register the product

    ## CHECK this code for correctness ...
    ## More than one documentation file is possible.
    ## Is it allowed to list more than one file for the 'help' key ?

    ## Bugzilla 18843 hack.
    ## When installing Tcl/DevKit also set keys for ActiveTcl,
    ## use a fixed version number for now. Later we will
    ## have to change this so that the build system is able
    ## to provide more information ... Maybe the build system
    ## should provide a list of keys and values, the latter
    ## having placeholders for install time information.

    ## Another hack: For Tcl/DevKit we know that there are
    ## two help files, and the second one is for ActiveTcl.

    ## Note: The upgrade kit does not have to set the ActiveTcl
    ##       keys, as it assumes that ActiveTcl is installed.


    set docfiles [list]
    foreach d $::AT(HELP) {
	if {[string equal $::AT(MODE) normal]} {
	    set d [file rootname $d]$::AT(tcl_VERSION).chm
	}

	lappend docfiles "$ndir\\doc\\$d"
    }

    set key "[soft_root]\\$::AT(Company)\\$::AT(PGROUP)"
    if {[catch {
	set_key  $key "CurrentVersion" "$::AT(VERSION)"
	set_key "$key\\$::AT(VERSION)" {} $ndir

	if {[llength $docfiles] > 0} {
	    set_key "$key\\$::AT(VERSION)\\Help" {} \
		    [lindex $docfiles 0]
	}
    } err]} { log "ERROR: $err" error }

    return $docfiles
}

proc set_komodo_search {installdir ndir bindir} {
    ## Create the standard information used by Komodo to look for interpreters.
    ## Not required for lite, we expect an existing ActiveTcl install
    ## and thus assume the existence of the association

    if {[string equal $::AT(MODE) pro]}        {return}

    global INSTALL_MODE ; # "admin" or "user"

    if {[string equal $INSTALL_MODE user]} {
	# Personal installation without admin privileges ...
	# (If this is a personal installation for a user with
	# admin privileges we go forward with our action.

	log " "
	log "\tFor a personal installation we assume that extending \"App Paths\" with" warning
	log "\ta reference to tclsh is not possible as this action requires access to"  warning
	log "\tprivileged keys in the registry."                                        warning
	log " "                                                                         warning
	log "\tConsequence: ActiveState Komodo may not find the tclsh/wish of this"     warning
	log "\t             installation."                                              warning
	log " "

	return
    }

    if {$::AT(DEBUG)} {
	set key {HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\tclshg.exe}
    } else {
	set key {HKEY_LOCAL_MACHINE\SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\tclsh.exe}
    }
    log "\tExtend \"App Paths\" with reference to tclsh"

    if {[catch {
	set_key $key {}   "$bindir\\tclsh$::AT(ND_VERG).exe"
	set_key $key Path "$bindir"
    } err]} { log "ERROR: $err" error }
    return
}

proc set_path_search {installdir ndir bindir} {
    ## Create the path information
    ## Not required for TDK, we expect an existing
    ## (Active)Tcl install and thus assume the existence of the
    ## association

    ## Untrue for TDK. For TDK we need the installdir/bin inteh path
    ## for cmdline access to the tools.

    global INSTALL_MODE ; # "admin" or "user"
    global systemEnv userEnv

    if {[string equal "Windows 95" $::tcl_platform(os)]} {
	# Windows 95/98/ME needs path added to AUTOEXEC.BAT
	log "NOTE: If you wish to have tclsh and wish in your"
	log "path, you must add \"$bindir\" to the PATH environment variable"
	log "in C:\\AUTOEXEC.BAT (takes effect upon restart)."
    } else {
	if {[catch {
	    log "\tAdding \"$bindir\" to your executable path ..."

	    # Cases:
	    # 1. Admin installation ...
	    #    a) bindir is in common PATH => Nothing to do
	    #    Extend common PATH.
	    # 2. Personal install ...
	    #    a) bindir is in common PATH   => Nothing to do.
	    #    b) bindir is in personal PATH => Nothing to do.
	    #       (That key may not exist !)
	    #    Extend personal PATH.

	    set curPath [r_get $systemEnv Path]

	    if {[regexp -nocase -- "***=$bindir;" $curPath]} {
		log "\t\t(\"$bindir\" already in your path)"
	    } else {
		switch -exact -- $INSTALL_MODE {
		    admin   {
			r_set $systemEnv Path "$bindir;$curPath" expand_sz
		    }
		    user {
			set curPath ""
			catch {set curPath [r_get $userEnv Path]}

			if {[regexp -nocase -- "***=$bindir;" $curPath]} {
			    log "\t\t(\"$bindir\" already in your path)"
			} else {
			    # Use %PATH% to enforce inheritance of
			    # existing path information from the
			    # global environment.

			    if {[string equal "" $curPath]} {
				set newPath "$bindir"
			    } else {
				set newPath "$bindir;$curPath"
			    }
			    if {![string match -nocase "*%PATH%*" $newPath]} {
				append newPath ";%PATH%"
			    }
			    r_set $userEnv Path $newPath expand_sz
			}
		    }
		    default {
			return -code error "PANIC.\
				Internal error of the installer.\
				This should not happen."
		    }
		}
	    }
	} err]} {
	    log "ERROR: $err" error
	}
    }

    if {![info exists ::INSTALL_PATHEXT] || $::INSTALL_PATHEXT} {
	if {[catch {
	    log "\tAdding \".tcl\" to your pathext ..."

	    # Cases:
	    # 1. Admin installation ...
	    #    a) .tcl is in common PATHEXT => Nothing to do
	    #    Extend common PATHEXT.
	    # 2. Personal install ...
	    #    a) .tcl is in common PATHEXT   => Nothing to do.
	    #    b) .tcl is in personal PATHEXT => Nothing to do.
	    #       (That key may not exist !)
	    #
	    #   Extend personal PATHEXT.
	    #
	    #   Bugzilla TclDevKit 20828: If it was empty, but the common
	    #   PATHEXT was not we have to copy the data from the common
	    #   PATHEXT to prevent the loss of these associations.

	    set nosystemdefault 0
	    if {[catch {
		set curSystemPath [r_get $systemEnv PATHEXT]
	    }]} {
		# The system environment variable PATHEXT was missing.
		# Fall back to hardwired defaults.
		set curSystemPath ".COM;.EXE;.BAT;.CMD"
		set nosystemdefault 1
	    }

	    if {[string match -nocase "*;.tcl*" $curSystemPath]} {
		log "\t\t(\".tcl\" already in your pathext)"
	    } else {
		switch -exact -- $INSTALL_MODE {
		    admin   {
			r_set $systemEnv PATHEXT "$curSystemPath;.tcl" sz
		    }
		    user {
			set curPersonalPath ""
			catch {set curPersonalPath [r_get $userEnv PATHEXT]}

			if {[string match -nocase "*;.tcl*" $curPersonalPath]} {
			    log "\t\t(\".tcl\" already in your pathext)"
			} else {
			    if {$curPersonalPath == ""} {
				# Bugzilla TclDevKit 20828
				#
				# Use common variable as source to prevent
				# the loss of the information contained
				# therein.

				set curPersonalPath "%PATHEXT%"
			    }
			    if {$nosystemdefault} {
				# Make sure that basic extensions are present if
				# the system PATHEXT was empty.
				append curPersonalPath ";" $curSystemPath
			    }
			    r_set $userEnv PATHEXT "$curPersonalPath;.tcl" expand_sz
			}
		    }
		    default {
			return -code error "PANIC.\
			    Internal error of the installer.\
			    This should not happen."
		    }
		}
	    }
	} err]} {
	    if {[string compare "Windows 95" $::tcl_platform(os)]} {
		# Only complain on Win2K/NT
		log "ERROR: $err" error
	    }
	}
    }

    catch {registry broadcast Environment}
    return
}

proc register_uninstaller {installdir ndir bindir} {
    ## Create the uninstall information

    global INSTALL_MODE

    set key "HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\$::AT(PGROUP) $::AT(MVERSION)"

    if {[string equal $INSTALL_MODE user]} {
	# Personal installation without admin privileges ...
	# (If this is a personal installation for a user with
	# admin privileges we go forward with our action.

	log " "
	log "\tFor a personal installation we assume that extending 'Settings | Control Panel | Add/Remove Programs'" warning
	log "\twith an uninstaller is not possible as this action requires access to"  warning
	log "\tprivileged keys in the registry."                                       warning
	log " "                                                                        warning
	log "\tThis means that the only way to uninstall this package is through the"  warning
	log "\tuninstall entry added to the start menu."                               warning
	log " "
	return
    }

    log "\tRegister uninstaller in system"

    if {[catch {
	set_key $key "DisplayName"    "$::AT(Company) $::AT(NAME) $::AT(MVERSION)"
	set_key $key "Publisher"      "$::AT(Company) Software Inc."
	set_key $key "DisplayVersion" "$::AT(MVERSION)"
	set_key $key "HelpLink"       $::AT(url,help)
	set_key $key "URLInfoAbout"   $::AT(url,about)
	set_key $key "URLUpdateInfo"  $::AT(url,update)

	set unapp [get_uninstallshell $ndir]

	set_key $key "UninstallString" "\"$unapp\""

	install_log [list REGKEY $key]
    } err]} {
	log "\tERROR: $err" error
    }
    return
}

proc get_uninstallshell {ndir} {
    # Windows specific code.
    return [file nativename \
		[file join \
		     [file dirname $::INSTALL_DIR/$::AT(LOG)] \
		     $::AT(UNINSTALLER)]]
}

proc teacup_tdk_setup {installdir} {
    # Works only partially if 'teacup' is not in the path.

    set thr [file join $::RUNTIME_DIR tap_help_repository]

    log " "
    log "TEA Package Management - TDK .tap Helper Repository."
    log " "
    log " Adding to TDK preferences ..."
    log " % $thr"

    tdk_extend_preferences pkgRepositoryList [list $thr]
    log " Ok"

    set cmd [teapot_client_tdk]
    if {$cmd eq ""} {
	log "No teacup in the PATH, not looking for installation repositories"
	return
    }

    log " "
    log "TEA Package Management - Locate shells and connected repositories."
    log " "

    set repos {}
    foreach p $::TCLPATHS {
	set shells [glob -nocomplain -directory $p tclsh*]
	if {![llength $shells]} continue
	foreach sh $shells {
	    if {[catch {
		set lines [split [eval [linsert [linsert $cmd 0 exec] end link info $sh]] \n]
	    } msg]} {
		log "  $msg" err
		continue
	    }
	    set first 1
	    foreach line $lines {
		if {![regexp -- {^Repository (.*)$} $line -> r]} continue
		lappend repos $r

		if {$first} {
		    set first 0
		    log " * [file nativename $sh]"
		}
		log "   = [file nativename $r]"
	    }
	}
    }

    if {![llength $repos]} return

    set repos [lrmdup $repos]

    log " "
    log " Adding to TDK preferences ..."
    foreach r $repos { log " % $r" }
    tdk_extend_preferences pkgRepositoryList $repos
    log " Ok"
    return
}

proc tdk_extend_preferences {key paths} {
    # Keep possibly existing paths

    set pathsenv [pref::prefGet $key GlobalDefault]

    foreach p $paths {lappend pathsenv $p}
    set paths [lrmdup $pathsenv]

    pref::prefSet   GlobalDefault $key $paths
    pref::groupSave GlobalDefault

    return
}


proc teacup_setup {installdir srcdir runtimedir} {
    ### ### ### ### ### #

    # No automatic teapot integration for debug builds. They do not contain the
    # teacup application used for this here.

    if {$::AT(DEBUG)} return

    ## TEAcup for a AT distributions, whatever Tcl version, also no
    ## differences between regular AT and ATEE.

    set url  http://teapot.activestate.com
    set tp  [teapot_client $installdir]
    set sh  [teapot_tclsh  $installdir]

    log " "
    log "TEA Package Management - Standard Repositories."
    log " "

    log "  TEAcup          @ $tp"

    # Notes: The first invokation of teacup in the series below will
    # automatically migrate (copy!) a shared configuration found under
    # "~/.teapot/config" to the new installation&platform-dependent
    # location under "~/.teapot/<INSTALL_CODE>/config".
    #
    # That behaviour is actually sensible.
    #
    # Of the 5 copied pieces of data (timeout, proxy, archive list,
    # default install repository location, index cache location) the
    # first three are something we really want to replicate in the new
    # install, as they are how the user wanted them and the new
    # installation will likely be setup in the same way.
    #
    # Copying the other/last two doesn't make sense, as they tie into
    # the location of the new installation. Doing so doesn't harm us
    # however, as the 'teacup default' and 'teacup cache on' commands
    # coming later below will overwrite the copied values with the
    # correct ones.

    log* "  Create repository"
    if {[catch {
	set res [split [exec $tp create [file nativename $installdir/lib/teapot]] \n]
	#set res [split [exec $tp create] \n]
	regsub {^Repository @ } [lindex $res 0] {} dir
	set dir [file nativename $dir]
    } err]} {
	log ": $err" error
	log "Aborting TEApot setup" warning
	return
    }
    log " @ $dir"

    # Connect repository and shell. Generates a file 'teapot-link.txt'
    # in the directory "installdir/lib/tcl8.X". Record this file in
    # the installer log so that it will be picked up and removed by
    # the uninstaller.

    if {[canonpath $installdir] ne [canonpath $runtimedir]} {
	package require fileutil

	# With installdir != runtimedir the teacup cannot be used to
	# make the connection. The runtimedir doesn't exist, so there
	# is no repository to connect to, and the shell cannot be run
	# either at this point, with the binary patching causing it to
	# look into the wrong place. We have to do everything
	# 'manually', using 'file commands and such.

	set rsh  [teapot_tclsh    $runtimedir]
	set rdir [file nativename $runtimedir/lib/teapot]

	log "  Setup shell     @ $rsh"
	log "                <-> $rdir"
	log "  MANUALLY"

	# 1. installdir/lib/tcl8.X/teapot-link.txt, point to future location of the repository.

	fileutil::appendToFile \
	    [file join $installdir lib tcl$::AT(tcl_VERSION) teapot-link.txt] \
	    $rdir\n

	# 2. lib/teapot/teapot.local-ma.config, point to future location of the shell.
	#    NOTE: Because teacup can be used only from there, everywhere else the shell
	#          can't be executed.

	fileutil::writeFile \
	    [file join $dir teapot.local-ma.config] \
	    [list [platform::identify] [list unix [list $sh 1]]]\n

	log "  Ok."
    } else {
	# Regular setup.

	subenv {
	    # We are dealing with a patched shell here, the env
	    # variables are not necessary. Actually having them is a
	    # bug because they cause the installed and patched shell
	    # we interogate here to still look in the distribution for
	    # its 'info library' directory, causing teacup to set the
	    # repository link information into the wrong location.

	    # NOTE: On windows the variables do not exist, and are not
	    # NOTE: needed either.

	    unset -nocomplain ::env(TCL_LIBRARY)
	    unset -nocomplain ::env(TK_LIBRARY)
	    unset -nocomplain ::env(LD_LIBRARY_PATH)
	    unset -nocomplain ::env(DYLD_LIBRARY_PATH)
	    unset -nocomplain ::env(SHLIB_PATH)
	    unset -nocomplain ::env(LIBPATH)

	    #See where the tclsh/teacup look for the repo link information.
	    #set ildir [exec $sh << {puts [info library];exit}]
	    #log "ildir($sh) = $ildir"
	    #log "has teapot.txt = [file exists $ildir/teapot.txt]"

	    log "  Setup shell     @ $sh"
	    log "                <-> $dir"
	    if {[catch {
		exec $tp link make $dir $sh
	    } err]} {
		log "   $err" error
		log "Aborting TEApot setup" warning
		return
	    }
	    log "  Ok."

	    #log "has teapot-link.txt = [file exists $ildir/teapot-link.txt]"
	}
    }

    install_log [list FILECOPY \
		     [file join $installdir lib tcl$::AT(tcl_VERSION) teapot-link.txt]]

    log* "  Default local   @ $dir"
    if {[catch {exec $tp default $dir} err]} {
	log " Error: $err" error
    } else {
	log " Ok."
    }

    log* "  Default network @ $url"
    if {[catch {exec $tp archive add $url} err]} {
	if {[string match "*already known*" $err]} {
	    log* " "
	    log* [string trim $err] url
	    log " "
	} else {
	    log " Error: $err" error
	}
    } else {
	log " Ok."
    }

    # Standard setup: cache repository index databases locally for
    # speed 'cache on' forces the index cache to the hardwired default
    # location per 'repository::sys::cachedir'. Anything in the
    # configuration is overwritten by this.
    log* "  Caching ... on"
    if {[catch {exec $tp cache on} err]} {
	log " Error: $err" error
    } else {
	log " Ok."
    }

    # Standard setup: Clear cache to force download on first use,
    # creating any directories possibly missing at that time.
    log* "  Caching ... cleared"
    if {[catch {exec $tp cache clear} err]} {
	log " Error: $err" error
    } else {
	log " Ok."
    }

    log "  Done"
    log " "

    # ### ### ### ######### ######### #########
    ## Time to look for a set of packages in the distribution to load
    ## into the local repository. If present, use teacup to load them
    ## into the repository.

    set preload [file join $srcdir payload lib teapot]

    if {
	($::INSTALL_SEED eq "merge") &&
	[file exists $preload]
    } {
	log " Merging seed repository with existing repository"

	set del 0
	if {[file system $preload] ne "native"} {
	    # Copy to disk for access by teacup.
	    set d [tempfile ATSEEDREPO_]
	    log " Making seed repository accessible on disk"
	    log " at $d"
	    file delete $d
	    file copy -force $preload $d
	    set preload $d
	    set del 1
	}

	set pd [tempfile ATSEEDPKG_]
	log " Retrieving packages from seed repository"
	log " to $pd"
	file delete $pd
	file mkdir  $pd
	set here [pwd]
	cd $pd
	pipe [list $tp get --at $preload]
	cd $here

	if {$del} {
	    log " Removing seed repository from disk"
	    log " Was at $preload"
	    file delete -force $preload
	}

	log " Adding seed packages to repository"
	log " at $dir"

	foreach packagefile [lsort -dict [glob -nocomplain -directory $pd *]] {
	    log "    Uploading $packagefile"
	    pipe [list $tp install --at $dir --force --no-follow $packagefile]
	}

	log " Removing seed packages from disk"
	log " Were at $pd"
	file delete -force $pd
    }

    return
}

proc pipe {cmd} {
    set pipe [open |$cmd r+]
    while {![eof $pipe]} {
	log "\t[gets $pipe]"
    }
    close $pipe
    return
}

proc tempfile {prefix} {
    package require fileutil
    return [fileutil::tempfile $prefix]
}

proc teapot_client_tdk {} {
    return [auto_execok teacup]
}

proc teapot_client {installdir} {
    set tp [file join $installdir bin teacup]
    if {[string equal "windows" $::tcl_platform(platform)]} {
	append tp .exe
    }
    return [file nativename $tp]
}

proc teapot_tclsh {installdir} {
    set sh [file join $installdir bin tclsh]
    if {1||[package vcompare $::AT(tcl_VERSION) 8.5] >= 0} {
	set tv $::AT(tcl_VERSION)
	if {[string equal "windows" $::tcl_platform(platform)]} {
	    set tv [string map {. {}} $tv]
	}
	append sh $tv
    }
    if {[string equal "windows" $::tcl_platform(platform)]} {
	append sh .exe
    }
    return [file nativename $sh]
}


proc patch_registry {installdir {group @@}} {
    ## patch registry and create shortcuts
    ## Extending the registry

    package require registry

    set ndir [file nativename $installdir]
    set bindir $ndir\\bin

    ## Create the file association
    ## Register the product
    ## Create the standard information used by Komodo to look for interpreters.
    ## Create the path information
    ## Create the uninstall information

    set_tcl_association  $installdir $ndir $bindir

    if {$::AT(MODE) ne "normal"} {
	# Pro distribution has to set associations for the
	# various TclDevKit project files. The icon is taken
	# from the basekit for Tk. As we have no regular wish.

	set launcher [list tcllauncher.exe]
	set xref     [list tclxref.exe]
	set icosrc   "$bindir\\tcllauncher.exe,0" \

	set_association \
	    TclDevKitProject "Tcl Dev Kit Project File" \
	    $launcher {} {.tpj .tap .tdk} $icosrc

	set_association \
	    TclDevKitXRefDB "Tcl Dev Kit Cross-Reference Database" \
	    $xref -gui {.txr} $icosrc
    }

    set docfiles [set_product $installdir $ndir $bindir]
    set_komodo_search    $installdir $ndir $bindir
    set_path_search      $installdir $ndir $bindir
    register_uninstaller $installdir $ndir $bindir

    ## Creating the shortcuts
    ## This action depends on the distribution the installer is run from.
    ##
    ## Differences: Availability of README's, Demos.

    if {$group eq "@@"} {
	set group "$::AT(Company) $::AT(dPGROUP) $::AT(MVERSION)"
    }

    ## FUTURE XXX FIXME : Declare the key in the 'data' files, use
    ##                    placeholders for the installation and demo
    ##                    directories.

    set smpdir [determine_startmenuprg]

    log_debug "Directory for start menu =\n\t$smpdir"

    set shortcut_items [list]
    # See [declare_menuentry] and [create_menu].

    if {[llength $docfiles] == 1} {
	set doc [lindex $docfiles 0]
	declare_menuentry "$::AT(dPGROUP) Help" $doc
    } elseif {[llength $docfiles] > 1} {
	foreach doc $docfiles {
	    set stem [file rootname [file tail $doc]]
	    if {$stem eq "TclDevKit"} {
		set stem "Tcl Dev Kit"
	    } elseif {$stem eq "ActiveTclHelp"} {
		# Bugzilla 25212
		set stem "ActiveTcl"
	    }
	    declare_menuentry "$stem Help" $doc
	}
    }

    # ### ### ### ######### ######### #########
    ## Menu uninstaller

    # Bugzilla 19781/19704 ... General changes: Use the new commands
    # to declare and create menu entries.

    set imgdir [file dirname $::INSTALL_DIR/$::AT(LOG)]/install_images

    # Use the -shortname because that's more portable
    set ico [file attributes "$imgdir/uninstall.ico" -shortname]

    set unapp [get_uninstallshell $ndir]

    declare_menuentry "Uninstall $::AT(dPGROUP) $::AT(MVERSION)" $unapp \
	-iconpath $ico

    ##
    # ### ### ### ######### ######### #########

    # ### ### ### ######### ######### #########
    ## Menu license

    # References to the licenses in the CHM files.
    # mk:@MSITStore:C:\Tcl\20tdk\doc\TclDevKit.chm::/Licensing.html
    # mk:@MSITStore:C:\Tcl\20tdk\doc\ActiveTclHelp.chm::/ActiveTcl8.3.4.3-html/at.license.html
    # Assume that they are not available anywhere.

    set doc [file nativename [file join $::INSTALL_DIR $::LICENSE]]
    # $::AT(NAME) - Bug 69693
    declare_menuentry "License" [get_editor_app] -args \"$doc\"
    ##
    # ### ### ### ######### ######### #########

    # ### ### ### ######### ######### #########
    ## Menu entries for ActiveTcl (shells, tkcon, readme, ...)

    if {$::AT(MODE) eq "normal"} {
	# Ensuring that the tcl/tk apps like tkcon and the demos use
	# the wish coming with the distribution, and not the wish of a
	# distribution installed afterward. IOW we bypass the .tcl
	# association.

	set instwish "$bindir\\$::wish_exe"

	declare_menuentry Wish$::AT(ND_VER) $instwish

	declare_menuentry Tclsh$::AT(ND_VER) \
	    "$bindir\\tclsh$::AT(ND_VERG).exe"

	if {1||($::AT(tcl_VERSION) eq "8.4")} {
	    declare_menuentry tkcon $instwish -args "\"$bindir\\tkcon.tcl\""
	}

	if {($::AT(tcl_VERSION) ne "8.4")} {
	    declare_menuentry {VFS Explorer} $instwish -args "\"$bindir\\tclvfse.tcl\""
	}

	if {0&&($::AT(tcl_VERSION) ne "8.4")} {
	    declare_menuentry guibuilder $instwish -args "\"$bindir\\guibuilder.tcl\""
	}

	declare_menuentry Readme "$ndir\\$::AT(README)"

	## Not yet. Requires PWD = demos/Snack to find its .wav
	## files. Installer does not handle such yet.
	## $instwish \"$demodir\\Snack\\widget.tcl\"" "Snack"

	# No demos yet: TclDOM, Tcllib, Tcl/DevKit, Soap, TclXML, TkHTML,
	#               Tktable
	#
	# Actually we have demos, a lot of demos for all of the
	# above. But either no central demo application, or the demos
	# have no GUI, making them pointless to use.

	set demodir [file nativename $::DEMO_DIR]
	set subtut Demos

	declare_menuentry "$subtut/Tk" $instwish \
	    -args "\"$demodir\\Tk$::AT(tcl_VERSION)\\widget.tcl\""

	if 0 {
	    # To reenable use different starkit, and ($::AT(tcl_VERSION) ne "8.4")
	    declare_menuentry "$subtut/Explore Starkit" \
		"$bindir\\tclvfse.tcl" -args "\"$bindir\\guibuilder.tcl\""
	}

	if {$::AT(tcl_VERSION) eq "8.4"} {
	    # Packages only in 8.4, thus their demos not in 8.5

	    declare_menuentry "$subtut/BWidget" $instwish \
		-args "\"$demodir\\BWidgets\\demo.tcl\""

	    declare_menuentry "$subtut/IWidgets" $instwish \
		-args "\"$demodir\\Iwidgets\\catalog.tcl\""

	    if {!$::AT(DEBUG)} {
		# Img demos not present in debug distributions.
		declare_menuentry "$subtut/Img" $instwish \
		    -args "\"$demodir\\Img\\demo.tcl\""
	    }

	    declare_menuentry "$subtut/TCom - Chart" $instwish \
		-args "\"$demodir\\Tcom\\chart.tcl\""
	}
    }


    if {$::AT(MODE) eq "pro"} {
	# Main TDK Apps
	declare_menuentry "TclApp"                "$bindir\\tclapp.exe"
	declare_menuentry "Package Editor"        "$bindir\\tclpe.exe"
	declare_menuentry "Compiler"              "$bindir\\tclcompiler.exe"

	if 0 {
	    # Avoid this entirely. XP SP3 has started to barf as well
	    # (compared to XP SP2). [Bug 77642].
	    if {!($::tcl_platform(os) eq "Windows NT"
		  && [package vsatisfies $::tcl_platform(osVersion) 6.0])} {
		# Avoid this on Vista as it barfs [Bug 76232]
		declare_menuentry "Checker" \
		    "mk:@MSITStore:$ndir\\doc\\TclDevKit.chm::/Checker.html"
	    }
	}

	declare_menuentry "Inspector"             "$bindir\\tclinspector.exe"
	declare_menuentry "Debugger"              "$bindir\\tcldebugger.exe"
	declare_menuentry "Cross Reference Tool"  "$bindir\\tclxref.exe"
	#declare_menuentry "VFS Explorer"          "$bindir\\tclvfse.exe"

	if {$::tcl_platform(os) ne "Windows 95"} {
	    # Install only for NT class machines
	    declare_menuentry "Tcl Service Manager" "$bindir\\tclsvc.exe"
	}

	declare_menuentry "Readme" "$ndir\\$::AT(README)"

	# TDK Demos
	set demodir [file nativename $::DEMO_DIR]
	set subtut Demos

	set pc "Profiling and Coverage"
	set sp "Debugging Subprocesses"
	set wr "Wrapping"
	set wa "$wr/Advanced"
	set wb "$wr/Bwidget Demo"
	set pe "Package Definitions"

	set instdbgr   "$bindir\\tcldebugger.exe"
	set insttclapp "$bindir\\tclapp.exe"

	declare_menuentry "$subtut/$pc/Hotspot Multiplication I" $instdbgr \
	    -args "\"$demodir\\TclDevKit\\TclDebugger\\Profiling\\mult_1.tpj\""

	declare_menuentry  "$subtut/$pc/Hotspot Multiplication II" $instdbgr \
	    -args "\"$demodir\\TclDevKit\\TclDebugger\\Profiling\\mult_2.tpj\""

	declare_menuentry "$subtut/$pc/Coverage Testsuite I" $instdbgr \
	    -args "\"$demodir\\TclDevKit\\TclDebugger\\Coverage\\csv.tpj\""

	declare_menuentry "$subtut/$pc/Coverage Testsuite II" $instdbgr \
	    -args "\"$demodir\\TclDevKit\\TclDebugger\\Coverage\\csv_new.tpj\""

	declare_menuentry $subtut/$sp $instdbgr \
	    -args "\"$demodir\\TclDevKit\\TclDebugger\\Spawning\\spawn.tpj\""

	declare_menuentry "$subtut/$wb/Create Starkit" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\BWidgets\\demo_starkit.tpj\""

	declare_menuentry "$subtut/$wb/Create Starpack" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\BWidgets\\demo_starpack.tpj\""

	declare_menuentry "$subtut/$wr/Img Demo -- As Starpack" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\Img\\demo.tpj\""

	declare_menuentry "$subtut/$wr/TkCon -- As Starkit" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\Tkcon\\demo.tpj\""

	declare_menuentry "$subtut/$wr/Klondike -- As Starpack" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\Klondike\\demo.tpj\""

	declare_menuentry "$subtut/$wa/TclDom -- Incorrect Wrap" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\Tcldom\\snoop-dom.tpj\""

	declare_menuentry "$subtut/$wa/TclDom -- Correct Complete Wrap" \
	    $insttclapp -args "-gui \"$demodir\\TclDevKit\\TclApp\\Tcldom\\snoop-dom-complete.tpj\""

	declare_menuentry "$subtut/$wr/Doctools Processor" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\TclApps\\dtp.tpj\""

	declare_menuentry "$subtut/$wr/Hexplode Game" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\TclApps\\hexplode.tpj\""

	declare_menuentry "$subtut/$wr/TkChat" $insttclapp \
	    -args "-gui \"$demodir\\TclDevKit\\TclApp\\TclApps\\tkchat.tpj\""

	if {$::tcl_platform(os) ne "Windows 95"} {
	    # Install only for NT class machines

	    declare_menuentry "$subtut/TclSvc Register ExprServer" \
		"$bindir\\tclsvc.exe" \
		-args "-script \"$demodir\\TclSvc\\exprserver.tcl\" \
			-servicename TclExpr -description Calculator"
	}

	# These items now rely on the presence of an ActiveTcl
	# installation providing the tclvfse application.
	declare_menuentry "$subtut/$wb/Explore Starkit" \
	    "tclvfse.tcl" -args "\"~\\bwdemo.tcl\""

	declare_menuentry "$subtut/$wb/Explore Starpack" \
	    "tclvfse.tcl" -args "\"~\\bwdemo.exe\""

	# Demos for TclXref and TclPE which do not depend on the
	# presence of an ActiveTcl installation.

	declare_menuentry "$subtut/Xref Encoders" "$bindir\\tclxref.exe" \
	    -args "-teapot \"$demodir\\TclDevKit\\TclXref\\encode\\teapot.txt\""

	declare_menuentry "$subtut/$pe/Encoders" "$bindir\\tclpe.exe" \
	    -args "-gui \"$demodir\\TclDevKit\\TclXref\\encode\\teapot.txt\""

	declare_menuentry "$subtut/$pe/Debugger Attach (Zip)" "$bindir\\tclpe.exe" \
	    -args "-gui \"$demodir\\TclDevKit\\TclPE\\package-tcldebugger_attach-1.4-tcl.zip\""

	declare_menuentry "$subtut/$pe/Debugger Attach (TM)" "$bindir\\tclpe.exe" \
	    -args "-gui \"$demodir\\TclDevKit\\TclPE\\package-tcldebugger_attach-1.4-tcl.tm\""
    }

    if {
	[string equal $::AT(MODE) normal] &&
	[string equal $::AT(tcl_VERSION) 8.4] &&
	[file exists $::DEMO_DIR/Expect/tkremotels.tcl]
    } {
	# Expect demos only for 8.4 installation, not 8.5, or where
	# Expect isn't available (like win64).

	log "\tExpect demos"

	set subtut Demos
	set demodir [file nativename $::DEMO_DIR]

	declare_menuentry "$subtut/Remote Ls" $instwish \
		-args "\"$demodir\\Expect\\tkremotels.tcl\""
    }

    # Now that the selection of menu entries is complete we can create them ...

    create_menu $group $shortcut_items $smpdir

    # ==================================================================

    ## Bugzilla item 18950
    ## Patch the *Config.sh files on Win* too

    if {[string equal $::AT(MODE) pro]}        {return}

    log "Patching configuration scripts ($::AT(tcl_VERSION))"

    set files  [glob -nocomplain -directory [file join $installdir lib] \
		    {t{cl,k}Config*.sh} \
		    "t{cl,k}$::AT(tcl_VERSION)/t{cl,k}Config*.sh"]

    ## log "pfx  = [list $::AT(BLD_PFXPATH)]"
    ## log "epfx = [list $::AT(BLD_EPFXPATH)]"

    # Changes in the insertvars causes u to get the paths as true
    # windows paths, i.e. with backward slashes. change here to,
    # create a variant with forward slashes now.

    set ekey [string map [list \\ /] $::AT(BLD_EPFXPATH)]
    set pkey [string map [list \\ /] $::AT(BLD_PFXPATH)]

    foreach f $files {
	log "*  [file nativename $f] ..."
	# Look and substitute both forms of a path, i.e.
	# with backward and forward slashes.

	patch_text_file $f $::AT(BLD_EPFXPATH) $installdir
	patch_text_file $f $ekey               $installdir
	patch_text_file $f $::AT(BLD_PFXPATH)  $installdir
	patch_text_file $f $pkey               $installdir
    }
    return
}


proc registry_root {} {
    global INSTALL_MODE ; # "admin" or "user"
    switch -exact -- $INSTALL_MODE {
	admin   {return HKEY_LOCAL_MACHINE}
	user    {return HKEY_CURRENT_USER}
	default {
	    return -code error "PANIC.\
		    Internal error of the installer.\
		    This should not happen."
	}
    }
}

proc soft_root {}  {return "[registry_root]\\SOFTWARE"}
proc class_root {} {return "[soft_root]\\Classes"}

proc determine_startmenuprg {} {
    # Query registry and environment to find the directory whose
    # contents are used to construct the start menu.

    # Bugzilla 19781 ... Admin/User Installation.

    global INSTALL_MODE ; # "admin" or "user"
    global env

    if {[string equal "Windows 95" $::tcl_platform(os)]} {
	# Win 95 this info cannot be found in HKLM, only HKCU
	set key "HKEY_CURRENT_USER\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\Explorer\\Shell Folders"
	set prg  "Programs"
	set sme  "Start Menu"
	set idx  USERPROFILE
    } else {
	set key "[soft_root]\\Microsoft\\Windows\\CurrentVersion\\Explorer\\Shell Folders"
	switch -exact -- $INSTALL_MODE {
	    admin {
		set prg "Common Programs"
		set sme "Common Start Menu"
		set idx  ALLUSERSPROFILE
	    }
	    user {
		set prg  "Programs"
		set sme  "Start Menu"
		set idx  USERPROFILE
	    }
	    default {
		return -code error "PANIC.\
			Internal error of the installer.\
			This should not happen."
	    }
	}
    }


    if {![catch {set path [r_get $key $prg]}]} {
	return $path
    }

    # Beware !! From here on the result is non-portable (i.e. not
    # internationalized) and will not work on non-US/non-English machines.

    log "Unable to find \"$key $prg\"" warning
    log "Falling back to non-i18n data" warning

    if {![catch {set path [r_get $key $sme]}]} {
	return [file join $path {Programs}]
    }

    log "Unable to find \"$key $sme\"" warning
    log "Falling back to environment variable \"$idx\"" warning

    if {[info exists env($idx)]} {
	return [file join $::env($idx) {Start Menu} Programs]
    }

    log "Unable to determine where the startmenu is located" warning
    return {}
}

# Bugzilla 19731/19704 ... Part of general reorganization of this part
# of the code. Creates the program group and start menu of the
# application. Distinguishes common and personal installation, and
# also if the installer was able to find the start menu in the file
# system.

proc create_menu {group shortcuts smpdir} {
    # shortcuts is a list of quartets defining
    # - (D) display name of shortcut
    # - (I) icon to use for shortcut
    # - (S) subfolder to place the shortcut in (= submenu)
    # - (C) command/file linked to by the shortcut
    # - (O) arguments to the command
    #
    # in this order. An empty icon signals to use the default one. An
    # empty directory signals that the entry goes into the main menu.
    #
    # 'list of quartets' means "shortcuts = {D I S C O D I S C O ...}".

    if {[catch {
	log "\tCreating program group \"$group\" ..."

	if {$smpdir == {}} {
	    # The system was unable to determine the location of the
	    # start menu in the file system.
	    return -code error \
		"unable to determine Start Menu Programs directory"
	}

	# We will require win32 1.1 (CreateShortcut) for this variant
	package require win32 1.1

	# The installer knows the location of the start menu in the
	# file system, and thus also the location of the directory
	# for any program group.

	set gdir [file join $smpdir $group]
	file mkdir $gdir
	install_log [list MKDIR $gdir]

	log_debug "Group directory = $gdir"

	foreach {name exe opts} $shortcuts {
	    # name may include subdir components
	    set lnkname [file join $gdir "$name.lnk"]
	    set subdir [file dirname $lnkname]
	    if {![file exists $subdir]} {
		# Make sure we capture all parents that don't exist
		# because mkdir implicitly creates the full structure,
		# and we need to delete anything we create [Bug #76274]
		set parent $subdir
		while {![file exists $parent]} {
		    install_log [list MKDIR $parent]
		    set parent [file dirname $parent]
		}
		file mkdir $subdir
	    }
	    log_debug "\tCreating link \"$lnkname\" to \"$exe\" ($opts) ..."
	    eval [list win32::CreateShortcut $lnkname $exe] $opts
	    install_log [list FILECOPY $lnkname]
	}
    } err]} {
	log "ERROR: $err" error
    }
    return
}

# Bugzilla 19731/19704 ... Part of general reorganization of this part
# of the code. Wrapper around data structures to make the declaration
# of new shortcuts easier to comprehend and maintain.

# was { name icon sub exe args }
proc declare_menuentry {name exe args} {
    upvar 1 shortcut_items si

    # Avoid space-based pathnames by using the shortname.
    if {0 && [string match "* *" $exe]} {
	set exe [file nativename [file attributes $exe -shortname]]
    }

    lappend si $name $exe $args
    return
}


proc do_finish {srcdir installdir runtimedir} {
    ## finish installation
    ##
    ## Windows: Registry ...
    ## Unix:    Patch tclsh, wish for installed location.

    log "Finishing the installation"

    if {[string equal $::AT(MODE) pro]} {
	# TDK installation. Look for (Active)Tcl installations and
	# make them available as standard TAP search paths. Also
	# standard local TEAPOT repositories.

	pristine_sentinel $installdir

	# ### ### ### ######### ######### #########
	##  -*- TDK -*- code sharing setup ---

	# On windows we have to mount the shared code to get access to
	# some packages we need for this (preferences). On unix the
	# installer is run by the starpack containing the common code
	# already, so no mounting is required.

	if {[string equal "windows" $::tcl_platform(platform)]} {
	    set src [file join $installdir bin tdkbase.mkf]
	    set dst [file join $installdir shared]

	    ::vfs::mk4::Mount $src $dst -readonly
	    lappend ::auto_path $dst/lib

	    # Need the platform specific directories as well.
	    foreach p [glob -nocomplain -directory $dst/lib P-*] {
		lappend ::auto_path $p
	    }

	    unset src
	    unset dst
	}

	locate_tcl_general
	get_activetcl_install
    }

    switch -glob -- $::tcl_platform(platform) {
	unix {
	    log "Patching the shells and libraries for the new location ..."
	    patch_shells $srcdir $installdir $runtimedir
	    patch_tut
	    patch_tools $installdir $runtimedir
	}
	win* {
	    log "Extending the registry ..."
	    patch_registry $installdir
	    patch_tut
	}
	default {
	    error "Unknown platform $tcl_platform(platform)"
	}
    }

    ### Repository integration
    ## AT has local installation repository, links to AS standard repository.

    if {[string equal $::AT(MODE) normal]} {
	teacup_setup $installdir $srcdir $runtimedir
    }
    if {[string equal $::AT(MODE) pro]} {
	teacup_tdk_setup $installdir
    }
    ### ### ### ### ### #

    log "Done"
    return
}

proc install_file {f dstdir} {
    file copy -force $f $dstdir
    install_log [list FILECOPY [file join $dstdir [file tail $f]]]
    return
}

proc do_install_modules {srcdir installdir demodir} {
    log "Installing $::AT(NAME) ..."

    install_all [file join $srcdir payload] $installdir

    foreach key {README LICENSE MANIFEST} {
	install_file [file join $srcdir $::AT($key)] $installdir
    }

    # Pro and normal only
    # Step II: Install the demos ...

    log "Installing demos ..."

    install_all [file join $srcdir pdemos] $demodir

    log "Done ..."
    return
}

proc contains_activetcl_installation {dir} {
    set ::AT(InstVersion) ""

    set found 0
    foreach pattern {
	README*
	bin/tclsh*
	bin/wish*
	bin/tkcon*
    } {
	if {[llength [eval glob -nocomplain -directory $dir $pattern]]} {
	    set found 1
	    break
	}
    }
    if {!$found} {
	# None of our regular files seem to be there.
	return 0
    }

    # Several key files were found.
    # Now check the contents of the README.

    # DANGER - The filename below is possibly bogus by now.
    set readme [file join $dir README.txt]
    if {[catch {open $readme} fh]} {
	# Couldn't open the README file...
	return 0
    }

    set readme [read $fh]
    close $fh
    # Return 0 if we didn't match, otherwise 1
    # Try a number of patterns, for various forms of a version number.

    foreach pattern {
	{ActiveTcl\S*\s*(\d+\.\d+\.\d+\.\d+)}
	{ActiveTcl\S*\s*(\d+\.\d+\.\d+)}
	{ActiveTcl\S*\s*(\d+\.\d+)}
    } {
	if {[regexp $pattern $readme -> ::AT(InstVersion)]} {
	    return 1
	}
    }
    return 0
}

proc contains_tdk_installation {dir} {
    set ::AT(InstVersion) ""

    set found 0
    foreach pattern {
	README-tdk*
	bin/tclapp*
	bin/tclpe*
    } {
	if {[llength [eval glob -nocomplain -directory $dir $pattern]]} {
	    set found 1
	    break
	}
    }
    if {!$found} {
	# None of our regular files seem to be there.
	return 0
    }

    # Several key files were found.
    # Now check the contents of the README.

    set readme [file join $dir README-tdk.txt]

    if {[catch {open $readme} fh]} {
	# Could not open the README for Tcl/DevKit ...
	return 0
    }

    set readme [read $fh]
    close $fh
    # Return 0 if we didn't match, otherwise 1


    foreach pattern {
	{TclDevKit\S*\s*(\d+\.\d+\.\d+\.\d+)}
	{TclDevKit\S*\s*(\d+\.\d+\.\d+)}
	{TclDevKit\S*\s*(\d+\.\d+)}
    } {
	if {[regexp $pattern $readme -> ::AT(InstVersion)]} {
	    return 1
	}
    }
    return 0
}


proc overinstall_allowed {} {
    # Are we allowed to overwrite the installation in the selected
    # installation directory?
    #
    # The logic is simple - if we satisfy the current version, which
    # means we are a >= version, then allow overinstallation

    # We are using AT(VERSION) because that is the version of the
    # payload to the installer. This is what has to be compatible with
    # the installed version. This assumes that InstVersion is for the
    # same thing as the payload is. This type of code will fail if
    # the payload is ActiveTcl, but InstVersion is for a TDK, for
    # example.

    return [package vsatisfies $::AT(VERSION) $::AT(InstVersion)]
}


# ----------------------------------------------
# INSTALLATION CHECK ROUTINES
# ----------------------------------------------

# It is high time to refactorize and place generic code into the
# library. This is all registry searching and then looking at the
# results. ...

# -- locate_component
#
# Check if the specified component/part has been registered,
# using the given installation mode. Version and key are returned
# through variables, and can be used for additional checks in/of
# the registry.
#
# Arguments
#    part       - Name of the component to look for.
#    admin      - Boolean flag. True indicates to look for an admin install.
#    versionvar - Variable to store found version info into.
#    keyvar     - Variable to store the used registry key into.
#
# Results
#    A boolean flag. True indicates that the component was found.
#
# Sideeffects
#    If the component was found the variables for version and key
#    are set. They are not touched otherwise.
#

proc locate_component {part admin versionvar keyvar} {
    # Assume existence of the component on unix
    if {![string equal "windows" $::tcl_platform(platform)]} {
	return -code error "This command is restricted to Windows"
    }

    package require registry

    if {$admin} {
	set root HKEY_LOCAL_MACHINE
    } else {
	set root HKEY_CURRENT_USER
    }

    set key "${root}\\SOFTWARE\\$::AT(Company)\\$part"

    if {[catch {registry get $key "CurrentVersion"} ver]} {
	return 0
    }

    upvar 1 $versionvar v $keyvar k
    set v $ver
    set k $key
    return 1
}

# --------------------------------------
# TDK installer, locate ActiveTcl installs.
# --------------------------------------

proc get_activetcl_install {} {
    # Locate a proper ActiveTcl installation we can use for the
    # demos (8.4.9.1+, has to provide the basekits, we also use AT
    # demos for ours).

    set ::ATdemodir    @atdemodir@    ; # Defaults which keep the placeholders if
    set ::ATinstalldir @atinstalldir@ ; # no suitable AT installation was found.

    # We assume that locate_tcl_general was already run.
    # Because it provides us with the information about
    # all Tcl installs found to whittle down for our use.
    # The global var ::TCLPATHS (bin paths to shells).

    set atp [locate_activetcl_installs]

    if {![llength $atp]} {
	log "The Tcl Dev Kit demos are impaired and will not work. They rely on the" warning
	log "presence of a suitable installation of ActiveTcl (8.4.9.1+). Such an"   warning
	log "installation was not found."                                            warning
	log " "

    } elseif {[llength $atp] > 1} {
	log "More than one suitable ActiveTcl installation was found for use by the" warning
	log "Tcl Dev Kit demos. Using the first found."                              warning
	log " "
    }
    if {[llength $atp]} {
	## See [xx] too for use.
	## DANGER. Not using the default location for demos during an
	## DANGER. ActiveTcl install will break the code below. To fix
	## DANGER. this an ActiveTcl installation has to record the
	## DANGER. location chosen for its demos in the registry, or
	## DANGER. somewhere else the TDK installer can find it. Maybe
	## DANGER. a small script under bin/ to get this information ?

	set atp            [file dirname [lindex $atp 0]]
	set ::ATdemodir    [file join $atp demos]
	set ::ATinstalldir $atp

	log "Patching the Tcl Dev Kit demos to use the ActiveTcl installation at"
	log "    [file nativename $atp]"
	log " "

	set b [file join $atp bin]

	log "Initialize TclApp Basekit search preferences"
	log "    [file nativename $b]"
	log " "

	pref::prefSet   GlobalDefault prefixPath  $b
	pref::prefSet   GlobalDefault interpPath  $b
	pref::prefSet   GlobalDefault iconPath    $b
	pref::groupSave GlobalDefault
    }
    return
}

proc locate_activetcl_installs {} {
    # Assume that locate_tcl_general was run before us.
    set paths {}

    log " "
    log "Locate an ActiveTcl 8.4.9.1+ installation for use by the Tcl Dev Kit demos."

    foreach p $::TCLPATHS {
	# Exclude installations without the basekits

	set kits [glob -nocomplain -directory $p base-t*]
	if {![llength $kits]} continue

	# Instead of looking at the filesystem we run the tclsh/wish
	# in question and ask it for its version, and if it truly is
	# an ActiveTcl shell.

	set shells [glob -nocomplain -directory $p tclsh*]
	if {![llength $shells]} {
	    set shells [glob -nocomplain -directory $p wish*]
	    if {![llength $shells]} continue
	}

	set atv [lindex [split [exec [lindex $shells 0] << {
	    if {[catch {package present ActiveTcl} msg]} {
		puts ""
	    } else {
		puts $msg
	    }
	    exit
	}] \n] end]

	# Not ActiveTcl
	if {$atv eq ""} continue

	# Not 8.4.9.1+
	if {[package vcompare $atv 8.4.9.1] < 0} continue

	lappend paths $p
	log "* [file nativename $p]"
    }

    log " "
    return $paths
}


# --------------------------------------
# TDK installer, locate Tcl installs.
# --------------------------------------

proc locate_tcl_general {} {

    init_preferences_projectinfo

    log " "
    log "Initialize TclApp TAP search preferences"

    set paths [locate_tcl_PATH]

    if {[string equal "windows" $::tcl_platform(platform)]} {
	foreach p [locate_tcl_registry] {
	    lappend paths $p
	}
    }

    # Look in the installation directory for a parallel installed AT.
    set paths [linsert $paths 0 [file join $::INSTALL_DIR bin]]

    # Paths = List of candidate paths to tcl shells. After 'check_tap'
    # it is a list of paths to viable library directories.

    set paths [lrmdup          $paths]

    set ::TCLPATHS $paths

    set paths [check_tap_paths $paths]
    set paths [lrmdup          $paths]

    if {![llength $paths]} {
	log "  Unable to initialize TAP search."                   warning
	log "  No Tcl installations were found."                   warning
	log "  Leaving possibly existing TAP search paths intact." warning
	return
    } else {
	foreach p $paths {log "\#  [file nativename $p]"}
    }

    # Now access the preferences of the just installed TDK, save our
    # paths into them, and at last save the preferences back. The
    # platform differences regarding the location of the preferences
    # are mostly handled by the preference core. But the base path is
    # something we have to construct here. For now. The 'projectInfo'
    # package might help in the future (See the init at the beginning
    # of the procedure).

    tdk_extend_preferences pkgSearchPathList $paths
    return
}

proc init_preferences_projectinfo {} {
    # We can do this, because we already have access to all internal
    # packages via the setup in 'do_finish', our (indirect) caller.
    package require projectInfo
    package require pref
    package require pref::devkit

    log " "
    log "Accessing TclDevKit preferences"
    log "    Root: $::projectInfo::prefsRoot"
    log "Location: $::projectInfo::prefsLocation"
    log " History: $::projectInfo::prefsLocationHistory"

    # Note: Now that we use a proper location history the preference
    # initialization will automatically pull existing settings
    # forward.

    log "Settings in the history are pulled forward."
    pref::devkit::init
    return
}

proc locate_tcl_PATH {} {
    global env

    # We do not care about duplicates. Our caller will remove any.

    if {[info exists env(PATH)]} {
	if {[string equal "windows" $::tcl_platform(platform)]} {
	    return [split $env(PATH) \;]

	} else {
	    return [split $env(PATH) :]
	}
    }

    return {}
}

proc locate_tcl_registry {} {
    if {![locate_component ActiveTcl [string equal $::INSTALL_MODE admin] v key]} {
	# No ActiveTcl found, therefore no paths to add either.
	return {}
    }

    # The key we get is the base key for _all_ versions of ActiveTcl
    # which are installed, not only the most current. Get all sub
    # keys, and then scan them for viable shells. Execute them as
    # usual to get the version and dir information.

    set res {}
    foreach k [registry keys $key *] {
	if {![catch {
	    set installdir [registry get $key\\$k ""]
	}]} {
	    lappend res [file join $installdir bin]
	}
    }

    return $res
}

proc check_tap_paths {paths} {
    global env
    set res {}

    set win [string equal "windows" $::tcl_platform(platform)]
    if {$win} {
	set patterns {
	    tclsh.exe  tclsh8*.exe
	    tclshg.exe tclsh8*g.exe
	    wish.exe   wish8*.exe
	    wishg.exe  wish8*g.exe
	}
    } else {
	set patterns {
	    tclsh  tclsh8.*
	    tclshg tclsh8.*g
	    wish   wish8.*
	    wishg  wish8.*g
	}
    }

    set next 0
    foreach p $paths {
	# Windows: Lexically normalize to catch more duplicates.
	if {$win} {set p [string tolower $p]}

	# Sort out missing paths
	if {![file exists $p]} continue

	foreach pt $patterns {
	    # Sort out paths which do not contain a proper tclsh/wish
	    set matches [glob -nocomplain -directory $p $pt]
	    if {![llength $matches]} continue

	    # Found proper shells. Run them to determine their library
	    # directory. Also sorts out all non-executables.

	    foreach shell $matches {
		set code [catch {exec $shell << {
		    puts "AS-TDK-V=[info tclversion]"
		    puts "AS-TDK-L=[info library]"
		    exit
		}} eres]
		if {$code == 0} {
		    # Bug 66029 - Use markers (s.a) and search for
		    # them, so that output generated by the user's
		    # .tclshrc will not interfere with us.
		    set v    {}
		    set tlib {}
		    foreach line [split $eres \n] {
			set line [string trim $line]
			if {[string match "AS-TDK-V=*" $line]} {
			    regexp {^AS-TDK-V=(.*)$} $line -> v
			}
			if {[string match "AS-TDK-L=*" $line]} {
			    regexp {^AS-TDK-L=(.*)$} $line -> tlib
			}
		    }
		    if {($v ne "") && ([package vcompare $v 8.4] >= 0)} {
			# Valid search directories.
			# (Both Tcl's lib dir, and lib above).
			lappend res $tlib
			lappend res [file dirname $tlib]

			# No need to check for other shells in the
			# current path.
			set next 1

			# No need to check other matches of the
			# same pattern.
			break
		    }
		} ;# else ignore problematic shells.
	    }
	    if {$next} {set next 0 ; break}
	}
    }

    # We do not care about duplicates. Our caller will remove any.
    return $res
}

proc lrmdup {list} {

    set win [string equal "windows" $::tcl_platform(platform)]

    array set tmp {}
    set res {}
    foreach p $list {
	set px $p
	if {$win} {
	    set px [file nativename $px]
	}
	if {[info exists tmp($px)]} continue
	lappend res $p
	set tmp($px) .
    }
    return $res
}

proc mmvers {version} {
    set lversion [split $version .]
    if {[llength $lversion] < 3} {return $version}
    return [join [lrange $lversion 0 1] .]
}

proc the_sentinel {installdir} {
    set sentinel [file join $installdir bin tdkbase]
    if {[string equal "windows" $::tcl_platform(platform)]} {
	append sentinel .mkf
    }
    return $sentinel
}

proc pristine_sentinel {installdir} {
    set sentinel [the_sentinel $installdir]
    # Take sentinel's current mtime and blank the low byte.
    # See License_V8/Runtime/License.h
    file mtime $sentinel [expr {[file mtime $sentinel] & 0xffffff00}]
    return
}

# --------------------------------------
# Process cmdline options for an automatic install without user
# intervention.
#
# Known options
#
# --directory DIR
#
#   Specify installation directory. Runtime directory is the
#   same. Demo directory is derived from that (DIR/demos).
#
#   ** Required ** This activates the automatic mode.
#
# --demo-directory DIR
#
#   Specify directory for demos. Optional. Ignored if
#   automatic mode was not activated with --directory.
#
# --runtime-directory DIR
#
#   Specify directory to patch in. Optional. Ignored if
#   automatic mode was not activated with --directory.
#   Ignored on Windows.
#
# --user-mode
#
#   Windows only, ignored on Unix. Install for the current
#   user. Default is to install for all users, i.e. admin-mode.
#

#puts *\t[join $pargv \n*\t]

if {[string equal "windows" $::tcl_platform(platform)]} {
    proc canonpath {p} {
	return [string tolower [file nativename $p]]
    }
} else {
    proc canonpath {p} {
	return $p
    }
}

proc shift {var} {
    upvar 1 $var v
    set res [lindex $v 0]
    set v [lrange $v 1 end]
    return $res
}

proc docmdline {} {
    global auto autocfg pargv

    set auto 0
    array set autocfg {}

    while {[llength $pargv]} {
	set o [shift pargv]
	if {![string match -* $o]} {
	    cmdlineerror $o
	}
	switch -exact -- $o {
	    --directory {
		set autocfg(installdir) [file nativename [file normalize [shift pargv]]]
		set auto 1
	    }
	    --demo-directory {
		set autocfg(demodir) [file nativename [file normalize [shift pargv]]]
	    }
	    --runtime-directory {
		set autocfg(runtimedir) [file nativename [file normalize [shift pargv]]]
	    }
	    --user-mode {
		set autocfg(mode) user
	    }
	    --keep-existing {
		set autocfg(repo) existing
	    }
	    --keep-seed {
		set autocfg(repo) seed
	    }
	    --merge-seed {
		set autocfg(repo)  merge
	    }
	    default {
		cmdlineerror $o
	    }
	}
    }
    if {$auto} {
	# We have to know the location of the license file, despite it
	# not being shown by the installer. For the start menu entry
	# refering to it.

	license_file

	if {![info exists autocfg(demodir)]} {
	    set autocfg(demodir) \
		[file join $autocfg(installdir) demos]
	}
	if {
	    ![info exists autocfg(runtimedir)] ||
	    [string equal $::tcl_platform(platform) windows]
	} {
	    # Runtime defaults to installdir, and on windows it has to be
	    # the installdir, whatever the user wants.
	    set autocfg(runtimedir) $autocfg(installdir)
	}
	if {![info exists autocfg(mode)]} {
	    set autocfg(mode) admin
	}
	if {![info exists autocfg(repo)]} {
	    set autocfg(repo) existing
	}
    }
}

# --------------------------------------

proc saveenv {} {
    global            env
    return [array get env]
}

proc restoreenv {h} {
    global      env
    array unset env *
    array set   env $h
    return
}

proc subenv {script} {
    set h [saveenv]
    set code [catch {uplevel 1 $script} r]
    restoreenv $h
    return -code $code $r
}

# --------------------------------------
