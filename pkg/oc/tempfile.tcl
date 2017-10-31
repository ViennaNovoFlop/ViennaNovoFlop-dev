# FILE: tempfile.tcl
#
# Instances of this class represent temporary files -- files created by
# an application for some purpose, but which should not remain after the
# application terminates.
#
# Last modified on: $Date: 2009-04-24 17:46:48 $
# Last modified by: $Author: dgp $

Oc_Class Oc_TempFile {

    common unique 0

    private proc DefaultDirectory {} {
        global tcl_platform env
        set result "."       ;# In case all else fails
        if {[string match windows $tcl_platform(platform)]} {
            if {[info exists env(TEMP)] && [file isdirectory $env(TEMP)] \
                    && [file writable $env(TEMP)]} {
                return $env(TEMP)
            }
            if {[info exists env(TMP)] && [file isdirectory $env(TMP)] \
                    && [file writable $env(TMP)]} {
                return $env(TMP)
            } 
            if {[info exists env(TMPDIR)] && [file isdirectory $env(TMPDIR)] \
                    && [file writable $env(TMPDIR)]} {
                return $env(TMPDIR)
            }
	    if {[info exists env(HOMEDRIVE)] && \
		    [info exists env(HOMEPATH)]} {
               set bdir "$env(HOMEDRIVE)$env(HOMEPATH)"
            } elseif {[info exists env(USERPROFILE)]} {
               set bdir "$env(USERPROFILE)"
            } else {
               catch {unset bdir}  ;# Safety
            }
	    if {[info exists bdir]} {
		set tdir [file join $bdir "Local Settings" TEMP]
		if {[file isdirectory $tdir] && [file writable $tdir]} {
		    return $tdir
		}
		set tdir [file join $bdir TEMP]
		if {[file isdirectory $tdir] && [file writable $tdir]} {
		    return $tdir
		}
		set tdir [file join $bdir TMP]
		if {[file isdirectory $tdir] && [file writable $tdir]} {
		    return $tdir
		}
                set tdir [file join $bdir AppData Local Temp]  ;# Windows Vista
		if {[file isdirectory $tdir] && [file writable $tdir]} {
		    return $tdir
		}
		if {[file isdirectory $bdir] && [file writable $bdir]} {
		    return $bdir
		}
	    }
            if {[file isdirectory C:/TEMP] && [file writable C:/TEMP]} {
                return C:/TEMP
            }
            if {[file isdirectory C:/TMP] && [file writable C:/TMP]} {
                return C:/TMP
            }

            if {[info exists env(PUBLIC)]} {  ;# Windows Vista
               set tdir $env(PUBLIC)
            } else {
               set tdir "C:/Users/Public"
            }
            set tdir [file join $tdir Documents]
            if {[file isdirectory $tdir] && [file writable $tdir]} {
               return $tdir
            }

            if {[file isdirectory C:/] && [file writable C:/]} {
                return C:/
            }
        } else { ;# Assume Unix
            if {[info exists env(TMP)] && [file isdirectory $env(TMP)] \
                    && [file writable $env(TMP)]} {
                return $env(TMP)
            } 
            if {[info exists env(TMPDIR)] && [file isdirectory $env(TMPDIR)] \
                    && [file writable $env(TMPDIR)]} {
                return $env(TMPDIR)
            }
            if {[info exists env(TEMP)] && [file isdirectory $env(TEMP)] \
                    && [file writable $env(TEMP)]} {
                return $env(TEMP)
            }
            if {[file isdirectory /tmp] && [file writable /tmp]} {
                return /tmp
            }
        }
        if {[file writable .]} {
            return .
        }
        error "Can't find a directory with permission\
		to write a temporary file"
    }

    const public variable stem = _
    const public variable extension = {}
    const public variable directory

    private variable absoluteName
    private variable fileHandle
    private variable cleanup

    Constructor {args} {
        eval $this Configure $args
        if {![info exists directory]} {
            if {[catch {[Oc_Config RunPlatform] GetValue \
                    path_directory_temporary} directory]} {
		set directory [$class DefaultDirectory]
	    }
        }

        set absoluteDir [file join [pwd] $directory]
        # Should force into absolute form.  Cheap check for now.
        if {![string match absolute [file pathtype $absoluteDir]]} {
	    # The Cygwin environment can show symptoms of multiple
	    # personality disorder.
	    global tcl_platform
	    if {![string match windows $tcl_platform(platform)] ||
		![string match absolute [file pathtype \
			     [file nativename $absoluteDir]]]} {
		error "Programming error:\
                           Temp filename $absoluteDir not absolute."
	    }
        }
        if {![file isdirectory $absoluteDir]} {
            set msg "Temporary directory '$absoluteDir' does not exist"
            error $msg $msg
        }
        if {![file writable $absoluteDir]} {
            set msg "Temporary directory '$absoluteDir' is not writable"
	    if {[lsearch -exact $args -directory] == -1} {
		# Caller did not name directory.  Our default is bad.
        	append msg "\nEdit the config file for your platform"
		append msg "\nto specify a writable temporary directory."
	    }
            error $msg $msg
        }
        set base [file join $absoluteDir $stem-[info hostname]-[pid]]
	foreach {fileHandle absoluteName unique} \
		[Oc_OpenUniqueFile -sfx $extension -start $unique \
                     -pfx $base] \
		break

	# Set up so that temp files get cleaned up on exit, unless
	# claimed by some caller.
	Oc_EventHandler New cleanup Oc_Main Exit [list $this Delete]
    }

    Destructor {
        if {[info exists fileHandle]} {
            catch {close $fileHandle}
        }
	if {[string compare Oc_Nop $cleanup]} {
	    $cleanup Delete
            if {[info exists absoluteName]} {
        	file delete $absoluteName
            }
	}
    }

    method Claim {} {
	$cleanup Delete
	set cleanup Oc_Nop
    }

    method AbsoluteName {} {
        return $absoluteName
    }

    method Handle {} {
        return $fileHandle
    }
}
