##########################################################################
# Class representing a scheduled event
##########################################################################
Oc_Class Oxs_Schedule {

    private array common index

    # Add an interface to change these
    private common events {Step Stage}

    # Is this 
    private array variable active
    private array variable frequency

    const public variable output
    const public variable destination

    Constructor {args} {
        eval $this Configure $args
        # Only one schedule per output, destination pair
        if {![catch {set index($output,$destination)} exists]} {
            $this Delete
            array set a $args
            set index($a(-output),$a(-destination)) $exists
            return $exists
        }
        set index($output,$destination) $this
        regsub -all {:+} $output : op
        regsub -all "\[ \r\t\n]+" $op _ op
        regsub -all {:+} $destination : d
        regsub -all "\[ \r\t\n]+" $d _ d
        trace variable active wu [subst { uplevel #0 {
            set "Schedule-$op-$d-active" \
            \[[list array get [$this GlobalName active]]]
        } ;# }]
        trace variable frequency wu [subst { uplevel #0 {
            set "Schedule-$op-$d-frequency" \
            \[[list array get [$this GlobalName frequency]]]
        } ;# }]
        foreach e $events {
            set active($e) 0
            set frequency($e) 1
        }

        # Remove the Schedule whenever output or destination goes away.
#       Oc_EventHandler New _ [Oxs_Output Lookup $output] Delete \
#               [list $this Delete] -groups [list $this]
        Oc_EventHandler New _ [Oxs_Destination Lookup $destination] Delete \
                [list $this Delete] -groups [list $this]
    }
    Destructor {
        Oc_EventHandler Generate $this Delete
        Oc_EventHandler DeleteGroup $this
        catch {unset index($output,$destination)}
    }
    proc Set {o d w e v} {
        $index($o,$d) Set$w $e $v
    }
    method Send {e} {
        upvar #0 [string tolower $e] count
        if {$count % $frequency($e)} {
            return
        }
        if {[catch {Oxs_Output Send $output $destination} msg]} {
            if {[catch {Oxs_Output Lookup $output}]} {
                $this Delete
            } else {
                return -code error -errorcode $errorCode $msg
            }
        }
    }
    method SetActive {event value} {
        set active($event) $value
        Oc_EventHandler DeleteGroup $this-$event
        if {$value} {
            Oc_EventHandler New _ Oxs $event [list $this Send $event] \
                    -groups [list $this-$event $this]
        }
    }
    method SetFrequency {event value} {
        set frequency($event) $value
    }
}

