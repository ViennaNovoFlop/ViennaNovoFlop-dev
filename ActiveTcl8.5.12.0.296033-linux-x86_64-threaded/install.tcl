# -*- tcl -*-
# Main installation script for ActiveTcl
# NO Tk
#
# Copyright 2001-2008 ActiveState Software Inc.
# All Rights Reserved.

package require Tcl 8.4

global argv argc
set pargv $argv ; set argv {} ; set argc 0

foreach d $auto_path {
    foreach a [glob -nocomplain -directory $d P-*] {
	lappend auto_path $a
    }
}

lappend auto_path [file join [file dirname [info script]] lib]

# Arguments: 1 - Location of the distribution.
# Asks the user for the location to install to.

proc main {} {
    # Default install directory
    if {$::auto} {
	# Auto accepts the license implicitly.

	set ::ACCEPT 1

	# Standard configuration when doing a fully automatic
	# installation. Partially influenced by command line
	# options.

	set ::INSTALL_DIR  $::autocfg(installdir)
	set ::DEMO_DIR     $::autocfg(demodir)
	set ::RUNTIME_DIR  $::autocfg(runtimedir)
	set ::INSTALL_MODE $::autocfg(mode)
	set ::INSTALL_SEED $::autocfg(repo)
    } else {
	set ::INSTALL_DIR [default_installdir]
    }

    # Default install error message
    set ::ERRMSG ""

    if {$::auto} {
	# Automatic configuration, skip the config pages, jump
	# directly to the execution of the installation.  We still
	# have to auto-click the 'NEXT' button to actually start the
	# install. For this we hook into 'wait_next'.

	set state 9
    } else {
	# Manual configuration, go through all pages.
	set state 1
    }

    while {1} {
	switch -exact $state {
	    1 { incr state [intro] }
	    2 { incr state [license] }
	    3 { incr state ; # we have no check_previous_install step here }
	    4 { incr state [get_installdir] }
	    5 { incr state [overinstall] }
	    6 { incr state [get_demodir] }
	    7 { incr state [get_runtimedir] }
	    8 { incr state [seed_repository] }
	    9 { incr state [install_ready] }
	    10 { exit 0 }
	    default {
		return -code error "Unknow run state \"$state\""
	    }
	}
    }
}

# ----------------------------------------------

proc intro {} {
    if {[string match "log *" $::WELCOME]} {
	eval $::WELCOME
    } else {
	puts stdout $::WELCOME
    }
    return [wait_next]
}

proc license {} {
    puts stdout [license_text]
    return [wait_nextc \
	    "Cancel         \[no\]  => \[RET]\nAccept License \[yes\] => 'A' >>" \
	    {cancel no} {accept yes}]
}

# ----------------------------------------------

proc get_installdir {} {
    if {[string length $::ERRMSG]} { puts stdout $::ERRMSG }
    puts stdout "Please specify the installation directory."

    wait_next "Path \[$::INSTALL_DIR\]: "
    if {[string length $::ANSWER]} {
	set ::INSTALL_DIR [file nativename [file normalize $::ANSWER]]
    }

    foreach {::INSTALL_DIR ::ERRMSG} [check_installdir $::INSTALL_DIR] {break}
    if {[string length $::ERRMSG]} {
	# Something wrong with dir, do this again
	return 0
    } elseif {[string length $::AT(InstVersion)]} {
	# That we are here means that overinstallation was/is allowed.
	# Moved to next step (overinstall check).

	set ::DEMO_DIR [default_demodir $::INSTALL_DIR]
	set ::RUNTIME_DIR $::INSTALL_DIR
	return 1
    } else {
	# Skip the overinstall check

	set ::DEMO_DIR [default_demodir $::INSTALL_DIR]
	set ::RUNTIME_DIR $::INSTALL_DIR
	return 2
    }
}

proc overinstall {} {
    puts stdout [overinstall_warning]
    return [wait_next_back]
}

# ----------------------------------------------

proc get_demodir {} {
    if {[string length $::ERRMSG]} { puts stdout $::ERRMSG }
    puts stdout "Please specify the directory for the demos."

    wait_next "Path \[$::DEMO_DIR\]: "
    if {[string length $::ANSWER]} {
	set ::DEMO_DIR [file nativename [file normalize $::ANSWER]]
    }
    foreach {::DEMO_DIR ::ERRMSG} [check_demodir $::DEMO_DIR] {break}
    if {[string length $::ERRMSG]} {
	# Something wrong with dir, do this again
	return 0
    } else {
	return 1
    }
}

# ----------------------------------------------

proc get_runtimedir {} {
    if {[string length $::ERRMSG]} { puts stdout $::ERRMSG }
    puts stdout "Please specify the runtime installation directory.\n
This is the directory the applications will see as their installation directory
when searching for packages and libraries, instead of the directory the files
were copied to. In most circumstances this is the same as the installation
directory chosen before."

    wait_next "Path \[$::RUNTIME_DIR\]: "
    if {[string length $::ANSWER]} {
	set ::RUNTIME_DIR [file nativename [file normalize $::ANSWER]]
    }
    foreach {::RUNTIME_DIR ::ERRMSG} [check_runtimedir $::RUNTIME_DIR] {break}
    if {[string length $::ERRMSG]} {
	# Something wrong with dir, do this again
	return 0
    } else {
	return 1
    }
}

# ----------------------------------------------

proc seed_repository {} {
    # state 8
    if {![file exists $::SCRIPT_DIR/payload/lib/teapot]} {
	# The distribution has no seed repository, ignore this panel
	# and skip to the next state. No seed repository will be
	# installed.
	set ::INSTALL_SEED existing
	return 1
    }
    if {![file exists $::INSTALL_DIR/lib/teapot]} {
	# The destination directory has no repository either, ignore
	# this panel, and skip to the next state. The seed repository
	# we have will be installed unconditionally, through simple
	# file copy (install_all handles that)
	set ::INSTALL_SEED existing
	return 1
    }

    if {![info exists ::INSTALL_SEED]} {
	set ::INSTALL_SEED existing
    }

    # The distribution has a seed repository, and the destination has
    # a repository as well. We have to ask the user about what to do.
    # Options given to her:

    # (1) Keep existing repository, ignore seed.
    # (2) Use seed, and remove existing repository to make way.
    # (3) Merge seed with existing repository (slow).

    puts "This distribution comes with a seed repository, which is found to be in"
    puts "conflict with an existing repository at the chosen destination for the"
    puts "installation. Please specify how to handle this conflict."

    while {1} {
	puts ""
	puts "Your options are"
	set n 0
	foreach text {
	    {Keep the existing repository, and ignore the seed repository}
	    {Replace the existing repository with the seed repository}
	    {Merge the contents of the seed repository into the existing repository}
	} code {
	    existing
	    seed
	    merge
	} {
	    if {$::INSTALL_SEED == $code} {
		puts "\[[incr n]\] (Default) $text"
	    } else {
		puts "\[[incr n]\]           $text"
	    }
	}
	puts -nonewline stdout "Please enter the number of your choice >> "
	flush stdout

	set res [gets stdin]
	if {$res != {}} {
	    switch -exact -- [string trimleft $res 0] {
		1 { set ::INSTALL_SEED existing ; break }
		2 { set ::INSTALL_SEED seed     ; break }
		3 { set ::INSTALL_SEED merge    ; break }
	    }
	} else {
	    # Keep intall_seed value, default chosen
	    break
	}
    }

    return 1
}

# ----------------------------------------------

proc install_ready {} {
    puts stdout "Press return to begin installation"
    puts stdout "     Installation Directory:  $::INSTALL_DIR"
    puts stdout "     Demos Directory:         $::DEMO_DIR"

    if {[string equal $::INSTALL_DIR $::RUNTIME_DIR]} {
	puts stdout "     Runtime Directory:       See Installation Directory"
    } else {
	puts stdout "     Runtime Directory:       $::RUNTIME_DIR"
    }

    set res [wait_next]
    if {$res < 0} {
	return $res
    }

    if {[file pathtype $::INSTALL_DIR] ne "absolute"} {
	set ::INSTALL_DIR [file join [pwd] $::INSTALL_DIR]
    }
    if {[file pathtype $::DEMO_DIR] ne "absolute"} {
	set ::DEMO_DIR [file join [pwd] $::DEMO_DIR]
    }
    if {[file pathtype $::RUNTIME_DIR] ne "absolute"} {
	set ::RUNTIME_DIR [file join [pwd] $::RUNTIME_DIR]
    }

    # Install all the files
    do_install_modules	$::SCRIPT_DIR $::INSTALL_DIR $::DEMO_DIR
    # Patch files
    do_finish		$::SCRIPT_DIR $::INSTALL_DIR $::RUNTIME_DIR

    set p [parting_message]
    if {$p != {}} {
	puts $p
    }
    teapot_message
    hype_message

    set ::noconfirm 1
    wait_next "Finish >>"
    exit 0
}

# ----------------------------------------------
# WAIT ROUTINES
# ----------------------------------------------

global noconfirm ; set noconfirm 0

proc abortretryignore {msg} {
    puts $msg

    while {1} {
	puts -nonewline "\nAbort, retry, or ignore \[a/R/i\] ?"
	flush stdout
	set answer [gets stdin]
	if {$answer == {}} {set answer c}
	switch -exact -- [string index $answer 0] {
	    a - A { return abort  }
	    r - R { return retry  }
	    i - I { return ignore }
	}
    }
    return

}

proc cancel {} {
    if {$::noconfirm} {exit 0}
    while {1} {
	puts -nonewline "\nDo you truly wish to cancel the installation \[y/N\] ? "
	flush stdout
	set answer [gets stdin]
	if {$answer == {}} {set answer n}
	switch -exact -- [string index $answer 0] {
	    y - Y {exit 0}
	    n - N {break}
	}
    }
    return
}

proc wait_next_back {} {
    while {1} {
	puts -nonewline stdout "\nCancel => C\nNext   => N\nBack   => B : "
	flush stdout

	set answer [gets stdin]
	switch -exact -- [string index $answer 0] {
	    C - C {cancel}
	    N - n {return  1}
	    B - b {return -1}
	}
    }
}

proc wait_next {{label {}}} {
    while {1} {
	if {$label == ""} {
	    set label "Cancel => C\nNext   => \[RET] >>"
	}

	puts -nonewline stdout "$label"
	flush stdout

	if {!$::auto} {
	    set res [gets stdin]
	} elseif {$::auto == 1} {
	    # State 9, install_ready, 1st wait_next.
	    # Start installation. Simulate kbd input.

	    set res \n ; puts -nonewline stdout $res ; flush stdout

	} elseif {$::auto == 2} {
	    # State 9, install_ready, 2nd wait_next.
	    # Finish installation. Simulate kbd input.

	    set res \n ; puts -nonewline stdout $res ; flush stdout
	}

	if {
	    $res != "" &&
	    [regexp Cancel $label] &&
	    [string match -nocase ${res}* "Cancel"]
	} {
	    cancel
	} else {
	    break
	}
    }
    puts stdout ""
    set ::ANSWER $res
    return 1
}

proc wait_nextc {{label {}} {cwords {}} {awords {}}} {
    # Like wait_next, but cancel is default !!

    if {$label == ""} {
	set label "Cancel => \[RET]\nNext   => N >>"
    }
    if {$cwords == {}} {set awords cancel}
    if {$awords == {}} {set awords next}

    while {1} {
	puts -nonewline stdout "$label"
	flush stdout

	set res [gets stdin]
	if {$res != {}} {
	    set found 0
	    foreach cw $cwords {
		if {[string match -nocase ${res}* $cw]} {set found 1 ; break}
	    }
	    if {$found} {cancel ; continue}

	    foreach aw $awords {
		if {[string match -nocase ${res}* $aw]} {set found 1 ; break}
	    }

	    if {$found} break
	} elseif {$res == "" && [regexp Cancel $label]} {
	    cancel
	}
    }
    puts stdout ""
    set ::ANSWER $res
    return 1
}

# ----------------------------------------------
# LOGGING ROUTINE
# ----------------------------------------------

proc log {msg {type ok}} {
    puts stdout "$msg"
}

proc log* {msg {type ok}} {
    puts -nonewline stdout "$msg"
}

proc cmdlineerror {text} {
    puts stderr "ERROR\tCommand line error"
    puts stderr "\tIllegal argument \"$text\""
    exit 1
}

# ----------------------------------------------
# GO TO IT
# ----------------------------------------------

set ::SCRIPT_DIR [file dirname [info script]]
set here [pwd] ; cd $::SCRIPT_DIR ; set ::SCRIPT_DIR [pwd] ; cd $here
source [file join $SCRIPT_DIR install_lib.tcl]
docmdline
main
