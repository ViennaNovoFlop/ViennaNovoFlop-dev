# FILE: account.tcl
#
# An object which plays the role of an account within an interpreter.
#
# Each host on the internet may have many accounts -- separate areas
# of storage each for the use of one or more users, usually with access
# controlled by some sort of security measure such as a password.  Each
# thread which runs on a host is owned by exactly one of the accounts on
# that host.  Thus, in a hierarchy of "ownership," accounts lie between
# hosts and threads.
#
# The Tcl event loop must be running in order for non-blocking socket
# communications to work.
#
# Last modified on: $Date: 2008-05-20 07:13:39 $
# Last modified by: $Author: donahue $
#

Oc_Class Net_Account {

    # Index of account objects by name -- used to avoid duplicate objects.
    array common index {}
    common dir

    ClassConstructor {
        set dir [file dirname [info script]]
    }

    # Name of host where the account this object represents resides
    const public variable hostname = localhost

    # Name of account this object represents
    const public variable accountname

    # Default initializer of the accountname.  Tries to determine
    # account name under which current thread is running.
    proc DefaultAccountName {} {
       global env tcl_platform
       set name {}
       if {[info exists env(OOMMF_ACCOUNT_NAME)]
           && [string length $env(OOMMF_ACCOUNT_NAME)]} {
          set name $env(OOMMF_ACCOUNT_NAME)
       } elseif {[info exists tcl_platform(user)]
           && [string length $tcl_platform(user)]} {
          set name $tcl_platform(user)
       } else {
          switch -- $tcl_platform(platform) {
             unix {
                if {[info exists env(USER)]
                    && [string length $env(USER)]} {
                   set name $env(USER)
                } elseif {[info exists env(LOGNAME)]
                          && [string length $env(LOGNAME)]} {
                   set name $env(LOGNAME)
                } elseif {![catch {exec whoami} check]
                          && [string length $check]} {
                   set name $check
                } elseif {![catch {exec id -u -n} check]
                          && [string length $check]} {
                   # whoami is deprecated by posix id command
                   set name $check
                } elseif {![catch {exec id} check]
                          && [string length $check]} {
                   # Solaris has a non-posix id command that
                   # doesn't support the -u and -n switches.
                   if {[regexp {uid=[0-9]+\(([^)]+)\)} $check dummy check2]
                       && [string length $check2]} {
                      set name $check2
                   }
                }
             }
             windows {
                # Check environment
                if {[info exists env(USERNAME)]
                    && [string length $env(USERNAME)]} {
                   set name $env(USERNAME)
                }
                if {![catch {file tail [Nb_GetUserName]} check]} {
                   set name $check
                }
                # Find some reliable entry in the registry?
             }
          }
       }
       if {[string match {} $name]} {
          # Ultimate fallback
          set name oommf
       }
       return $name
    }

    # Milliseconds to wait for account thread to start on localhost
    public variable startWait = 60000
    # The ID of the current timeout (returned from 'after' command)
    private variable timeoutevent

    # The host which on which this account resides
    private variable host

    # Our connection to the actual account thread
    private variable connection

    # Port on host on which account thread listens
    private variable port
    
    # Socket for receiving ping from localhost OOMMF account server
    private variable listensocket

    # Send id
    private variable id = 0

    # Reply message
    private variable reply = {}

    # Is object ready to relay messages to and from the host?
    private variable ready = 0

    private variable hostDeathHandler

    private variable protocolversion

    public variable gui = 0 {
	package require Tk
	wm geometry . {} ;# Make sure we have auto-sizing root window
        if {$gui} {   
            if {[info exists myWidget]} {
                return
            }
            Net_AccountGui New myWidget $parentWin -account $this
            pack [$myWidget Cget -winpath] -side left -fill both -expand 1
        } else {
            if {![info exists myWidget]} {
                return
            }
            $myWidget Delete
            unset myWidget
	    $parentWin configure -width 1 -height 1  ;# Auto-resize
	    $parentWin configure -width 0 -height 0
        }
        set btnvar $gui
    }    
    public variable parentWin = .
    private variable myWidget
    private variable btnvar

    method ToggleGui { pwin } {
        set parentWin $pwin
        $this Configure -gui $btnvar
    }

    # Constructor is a non-blocking operation.  After it returns an instance
    # id, the Ready method of that instance will return a boolean indicating
    # whether the object is ready.  If not, then when it becomes ready, it
    # will generate an event with id 'Ready'.  If the object cannot become
    # ready, it will delete itself, generating an event with is 'delete'.
    # Users of instances of this class should set up to handle those two
    # events, when the Ready method returns false.
    #
    Constructor { args } {
        eval $this Configure $args
        if {![info exists accountname]} {
            set accountname [$class DefaultAccountName]
        }
        # Enforce only one object per account
        if {![catch {set index($hostname,$accountname)} existingaccount]} {
            unset hostname
            unset accountname
            $this Delete
            array set argarray $args
            catch {unset argarray(-hostname)}
            catch {unset argarray(-accountname)}
            eval $existingaccount Configure [array get argarray]
            return $existingaccount
        }
        array set index [list $hostname,$accountname $this]
	$this Restart
    }

    method Restart {} {
	set ready 0
        if {[catch {Net_Host New host -hostname $hostname} msg]} {
            catch {unset host}
            set msg "Can't construct host object.\n$msg"
            error $msg $msg
        } else {
            Oc_EventHandler New hostDeathHandler $host Delete \
		    [list $this FatalError "Lost connection to host\
		    $hostname"] -oneshot 1 -groups [list $this]
	    $this HandlePing
        }
    }

    private method HandlePing {args} {
        if {[llength $args]} {
            Oc_Log Log "$this: Received ping!" status $class
            close [lindex $args 0]
        }
        if {[$host Ready]} {
	    $this LookupAccountPort
        } else {
	    Oc_EventHandler New _ $host Ready [list $this LookupAccountPort] \
		    -oneshot 1 -groups [list $this]
        }
    }

    method LookupAccountPort {args} {
        set qid [$host Send lookup $accountname]
        Oc_EventHandler New _ $host Reply$qid \
                [list $this GetAccountPort $qid] \
                -groups [list $this $this-$qid]
        Oc_EventHandler New _ $host Timeout$qid \
                [list $this NoAnswerFromHost $qid] \
                -groups [list $this $this-$qid]
    }

    method GetAccountPort { qid } {
        Oc_EventHandler DeleteGroup $this-$qid
        set reply [$host Get]
        if {[llength $reply] < 2} {
            set msg "Bad reply from host $hostname: $reply"
            error $msg $msg
        }
        switch -- [lindex $reply 0] {
            0 {
                set port [lindex $reply 1]
                if {[info exists listensocket]} {
                    close $listensocket
                    unset listensocket
                }
                if {[info exists timeoutevent]} {
                    after cancel $timeoutevent
                    unset timeoutevent
                }
                $this Connect
            }
            default {
                # Host returns an error - assume that means lookup failed
                if {[string match localhost [$host Cget -hostname]]} {
                    if {![info exists listensocket]} {
                        $this StartLocalAccount
                    }
                } else {
                    $this FatalError "No account $accountname registered\
			    with host [$host Cget -hostname]"
                }
            }
        }
    }

    method NoAnswerFromHost { qid } {
        Oc_EventHandler DeleteGroup $this-$qid
        # Host not responding -- give up.
        $this FatalError "No reply from host [$host Cget -hostname]"
    }

    # Launch an OOMMF host server on localhost
    private method StartLocalAccount {} {
        set acctthread [file join $dir threads account.tcl]
        if {[file readable $acctthread]} {
            # Set up to receive ping from account thread or timeout
            set listensocket [socket -server [list $this HandlePing] \
		    -myaddr 127.0.0.1 0]  ;# Force localhost for now
            set listenport [lindex [fconfigure $listensocket -sockname] 2]
            set timeoutevent [after $startWait $this StartTimeout]
            Oc_Log Log "Starting OOMMF account server for $accountname\
		    on localhost" status $class
            # Ought to redirect std channels somewhere
	    if {[info exists lastoid]} {
		Oc_Application Exec {omfsh 1.1} \
		$acctthread -tk 0 0 $listenport $lastoid \
		> [[Oc_Config RunPlatform] GetValue path_device_null]  \
		2> [[Oc_Config RunPlatform] GetValue path_device_null] &
	    } else {
		Oc_Application Exec {omfsh 1.1} \
		$acctthread -tk 0 0 $listenport \
		> [[Oc_Config RunPlatform] GetValue path_device_null]  \
		2> [[Oc_Config RunPlatform] GetValue path_device_null] &
	    }
        } else {
            set msg "Install error: No OOMMF account server available to start"
            error $msg $msg
        }
    }

    method StartTimeout {} {
        $this FatalError "Timed out waiting for OOMMF account server\
		to start\nfor $accountname on localhost"
    }

    method Connect {} {
        if {[catch {
                Net_Connection New connection -hostname $hostname -port $port
                } msg]} {
            catch {unset connection}
            $this DeregisterBadAccount "Can't create connection:\n$msg"
        } else {
	    # Using [socket -async], it's possible that we haven't really
	    # yet established a valid connection.  So, until the connection
	    # is "Ready", handle all Net_Connection destruction as though
	    # caused by failure to connect.
            set h [Oc_EventHandler New _ $connection Delete \
		    [list $this DeregisterBadAccount \
		    "Connection to $hostname:$port failed"] \
                    -oneshot 1 -groups [list $this]]
            Oc_EventHandler New _ $connection Ready \
                    [list $this ConnectionReady $h] -oneshot 1 \
		    -groups [list $this]
        }
    }

    method ConnectionReady {handler} {
	$handler Delete
	Oc_EventHandler New _ $connection Delete [list $this Restart] \
		-oneshot 1 -groups [list $this]
        Oc_EventHandler New _ $connection Readable \
                [list $this VerifyOOMMFAccount] -oneshot 1 -groups [list $this]
    }

    method VerifyOOMMFAccount {} {
        if {![regexp {OOMMF account protocol ([0-9.]+)} [$connection Get] \
                match version]} {
            # Not an OOMMF account
            $this DeregisterBadAccount \
                    "No OOMMF account server answering on $hostname:$port"
            return
        }
        set protocolversion $version
	# Update when account server protocol version changes!
	if {![package vsatisfies $protocolversion 1]} {
	    # Old version of account server running!
            $this DeregisterBadAccount \
                    "Incompatible OOMMF account server is running\
		    at $hostname:$port."
            return
	}
        Oc_Log Log "OOMMF account $hostname:$accountname ready" status $class
        Oc_EventHandler New _ $connection Readable [list $this Hear] \
                -groups [list $this]
        set ready 1
	# Once we've connected to the account server, we don't need the
	# host server anymore.  Let it die and release its sockets.
        $hostDeathHandler Delete
        $host Delete
	unset host

	# Now get an OOMMF Id from the account server to represent our
	# application
	if {[info exists oid]} {
	    set qid [$connection Query \
		    newoid [Oc_Main GetAppName] [pid] [Oc_Main GetStart] $oid]
	} else {
	    set qid [$connection Query \
		    newoid [Oc_Main GetAppName] [pid] [Oc_Main GetStart]]
	}
        Oc_EventHandler New _ $connection Reply$qid [list $this Getoid $qid] \
                -oneshot 1 -groups [list $this $this-$qid]
        Oc_EventHandler New _ $connection Timeout$qid [list $this Nooid $qid] \
                -oneshot 1 -groups [list $this $this-$qid]
    }

    private variable oid
    private variable lastoid
    method OID {} {
	return $oid
    }

    method Getoid {qid} {
        Oc_EventHandler DeleteGroup $this-$qid
        set result [$connection Get] 
	if {[lindex $result 0]} {
	    $this FatalError "Did not receive OID: [lindex $result 1]"
	    return
	}
	set oid [lindex $result 1]
	if {![info exists lastoid]} {
	    set lastoid $oid
	    incr lastoid
	}
	if {[llength [info commands wm]]} {
	    #wm title . "<$oid> [Oc_Main GetAppName]|$lastoid"
	    wm title . "<$oid> [Oc_Main GetAppName]"
	    wm iconname . [wm title .]
	}
        Oc_EventHandler Generate $this Ready
    }

    method Nooid {qid} {
        Oc_EventHandler DeleteGroup $this-$qid
	$this FatalError "No response to OID request"
    }

    method DeregisterBadAccount {msg} {
        Oc_Log Log "Sending deregistration request for non-answering account" \
		status $class
        set qid [$host Send deregister $accountname $port]
        Oc_EventHandler New _ $host Reply$qid \
                [list $this DeregisterBadAccountReply $msg $id] \
                -groups [list $this $this-$qid]
        Oc_EventHandler New _ $host Timeout$qid \
                [list $this NoAnswerFromHost $qid] \
                -groups [list $this $this-$qid]
    }

    method DeregisterBadAccountReply {msg qid} {
        Oc_EventHandler DeleteGroup $this-$qid
        set reply [$host Get]
        switch -- [lindex $reply 0] {
            0 {
                if {[string match localhost [$host Cget -hostname]]} {
                    $this StartLocalAccount
                } else {
                    $this FatalError $msg
                }
            }
            default {
                $this FatalError "Unable to correct bad account server\
			registration after:\n$msg"
            }
        }
    }

    # Handle 'tell' messages from connection
    method Hear {} {
        set reply [$connection Get]
	if {[string compare [lindex $reply 0] private]} {
	    Oc_EventHandler Generate $this Readable
	    return
	}
	set lastoid [lindex $reply 1]
if {0 && [llength [info commands wm]]} {
    regexp {^([^|]*)[|]?.*$} [wm title .] -> p _
    wm title . "$p|$lastoid"
}
    }

    method Receive { rid } {
        Oc_EventHandler DeleteGroup $this-$rid
        set reply [$connection Get] 
        Oc_EventHandler Generate $this Reply$rid
    }

    method Get {} {
        return $reply
    }

    # Send a command to the host
    method Send { args } {
        if {!$ready} {
            set msg "account $hostname:$accountname not ready"
            error $msg $msg
        }
        incr id
	$this Resend $id $args
        return $id
    }

    private method Resend { x argList } {
        Oc_EventHandler DeleteGroup $this-$x
        set qid [eval $connection Query $argList]
        Oc_EventHandler New _ $connection Reply$qid [list $this Receive $x] \
                -oneshot 1 -groups [list $this $this-$x]
        Oc_EventHandler New _ $connection Timeout$qid [list $this Timeout $x] \
                -oneshot 1 -groups [list $this $this-$x]
        Oc_EventHandler New _ $this Ready [list $this Resend $x $argList] \
                -oneshot 1 -groups [list $this $this-$x]
    }

    # This routine silently cleans up after a query to the account server
    # times out with no response.
    method Timeout { rid } {
        Oc_EventHandler DeleteGroup $this-$rid
        Oc_EventHandler Generate $this Timeout$rid
    }

    method Ready {} {
        return $ready
    }

    method FatalError { msg } {
        global errorInfo
        set errorInfo "No stack available"
        Oc_Log Log "$this: $msg" warning $class
        $this Delete
    }
 
    Destructor {
        set ready 0
        Oc_EventHandler Generate $this Delete
        Oc_EventHandler DeleteGroup $this
	if {[info exists oid]} {
	    if {[llength [info commands wm]]} {
		wm title . [Oc_Main GetInstanceName]
		wm iconname . [wm title .]
	    }
	    set qid [$connection Query freeoid $oid]
	    Oc_EventHandler New _ $connection Reply$qid \
		    [list $connection Delete] -groups [list $connection]
	    Oc_EventHandler New _ $connection Timeout$qid \
		    [list $connection Delete] -groups [list $connection]
        } elseif {[info exists connection]} {
            $connection Delete
        }
        if {[info exists timeoutevent]} {
            after cancel $timeoutevent
        }
        if {[info exists listensocket]} {
            close $listensocket
        }
        if {[info exists myWidget]} {
            $this Configure -gui 0
        }
        if {[info exists hostname] && [info exists accountname]} {
            unset index($hostname,$accountname)
        }
    }

}
