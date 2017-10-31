# FILE: avf2ovf.tcl
#
# This file must be evaluated by mmDispSh

package require Oc 1.2.0.4        ;# Oc_TempName
package require Mmdispcmds 1.1.1.0
Oc_ForceStderrDefaultMessage
catch {wm withdraw .}

Oc_Main SetAppName avf2ovf
Oc_Main SetVersion 1.2.0.4

Oc_CommandLine Option console {} {}

proc FlipStringCheck { str } {
    set okaystr 0
    catch {
	if {[regexp {^[-+]?([xyz]):[-+]?([xyz]):[-+]?([xyz])$} \
		$str dummy c1 c2 c3] \
		&& ![string match $c1 $c2] \
		&& ![string match $c1 $c3] \
		&& ![string match $c2 $c3]} {
	    set okaystr 1  ;# Good flip string
	}
    }
    return $okaystr
}

proc ReadFile { filename } {
    if {[string match {} $filename]} {
        # Copy stdin to a temporary file
        set tempfile [Oc_TempName any2ovf]
        set outchan [open $tempfile "w"]
        fconfigure stdin -buffering full -translation binary
        fconfigure $outchan -buffering full -translation binary
        fcopy stdin $outchan
        close $outchan
        set workname $tempfile
    } else {
        if {![file exists $filename]} {
            return -code error "Input file \"$filename\" does not exist."
        }
        if {![file readable $filename]} {
            return -code error "Input file \"$filename\" is not readable."
        }

        # Apply filter (based on filename extension), if applicable.
        set tempfile [Nb_InputFilter FilterFile $filename decompress ovf]
        if { [string match {} $tempfile] } {
            # No filter match
            set workname $filename
        } else {
            # Filter match
            set workname $tempfile
        }
    }

    # Read file
    ChangeMesh $workname 0 0 0 -1
    
    # Delete tempfile, if any
    if { ![string match {} $tempfile] } {
	catch { file delete $tempfile }
    }

    if {[string match "Vf_EmptyMesh" [GetMeshType]]} {
	return -code error \
		"ERROR loading input file \"$filename\": Empty mesh"
    }

}

Oc_CommandLine Option flip {
    {flipstr { FlipStringCheck $flipstr } {is like z:-y:x (axes permutation)}}
} {
    upvar #0 flipstr globalflipstr
    regsub -all -- {[+]} $flipstr {} globalflipstr
} {Coordinate tranform}
set flipstr x:y:z

proc IsFloat { x } {
    if {[regexp -- {^-?[0-9]+\.?[0-9]*((e|E)(-|\+)?[0-9]+)?$} $x] || \
            [regexp -- {^-?[0-9]*\.?[0-9]+((e|E)(-|\+)?[0-9]+)?$} $x]} {
        if {![catch {expr $x}]} { return 1 }
    }
    return 0
}

proc ClipItemCheck { item } {
    # Format check
    if {[string match "-" $item] || [IsFloat $item]} {
	return 1
    }
    return 0
}

proc ClipStringCheck { clipstr } {
    set okaystr 0
    catch {
	if {[llength $clipstr]==0} {
	    set okaystr 1  ;# No clipping
	} elseif {[llength $clipstr]==6} {
	    set okaystr 1
	    # Format check
	    foreach item $clipstr {
		if {![ClipItemCheck $item]} {
		    set okaystr 0
		    break
		}
	    }
	    # Relative order check
	    for {set i 0} {$i<3} {incr i} {
		set min [lindex $clipstr $i]
		set max [lindex $clipstr [expr {$i+3}]]
		if {![string match "-" $min] && \
			![string match "-" $max] && \
			$min>$max} {
		    set okaystr 0
		    break
		}
	    }
	}
    }
    return $okaystr
}

Oc_CommandLine Option clip {
    {xmin { ClipItemCheck $xmin } {is minimum x}}
    {ymin { ClipItemCheck $ymin } {is minimum y}}
    {zmin { ClipItemCheck $zmin } {is minimum z}}
    {xmax { ClipItemCheck $xmax } {is maximum x}}
    {ymax { ClipItemCheck $ymax } {is maximum y}}
    {zmax { ClipItemCheck $zmax } {is maximum z}}
} {
    global clipstr
    set clipstr [list $xmin $ymin $zmin $xmax $ymax $zmax]
} "\n\tClipping box (optional)"
set clipstr {}

Oc_CommandLine Option keepbb {} {
    global cliprange; set cliprange 0
} {Don't clip bounding box.}
set cliprange 1

Oc_CommandLine Option format {
	{format {regexp {text|b4|b8} $format} {is one of {text,b4,b8}}}
    } {
	upvar #0 format globalformat; set globalformat $format
} {Format of output file: text, 4-, or 8-byte binary}
set format text

Oc_CommandLine Option grid {
	{type {regexp {^(reg|irreg)} $type} {is one of {reg,irreg}}}
    } {
       global grid ; regexp -- {^(reg|irreg)} $type grid
} {Grid type of output file: regular or irregular}
set grid reg

Oc_CommandLine Option info {} {global info; set info 1} {Print mesh info}
set info 0

Oc_CommandLine Option q {} {global quiet; set quiet 1} {Quiet}
set quiet 0

Oc_CommandLine Option [Oc_CommandLine Switch] {
    {{infile optional}  {} {Input vector field file (default: stdin)}}
    {{outfile optional} {} {Output vector field file (default: stdout)}}
} {
    global inputfile outputfile
    set inputfile $infile
    set outputfile $outfile
} {End of options; next argument is input file}
set inputfile {}
set outputfile {}

Oc_CommandLine Parse $argv

if {[string compare "-" $inputfile]==0} { set inputfile {} }
if {[string compare "-" $outputfile]==0} { set outputfile {} }

array set formatmap {
	text	text
	b4	binary4
	b8	binary8
}
array set gridmap {
	reg	regular
	irreg	irregular
}

# Check clip string is sensible
if {![ClipStringCheck $clipstr]} {
    puts stderr "ERROR: Bad clip request"
    exit 10
}

# Read file into mesh
if {![string match {} $inputfile]} {
    set inputfile [Oc_DirectPathname $inputfile]
}
if {[catch {ReadFile $inputfile} errmsg]} {
    puts stderr $errmsg
    exit 20
}

set title [GetMeshTitle]
set desc [GetMeshDescription]

if {$info} {
    if {[string match {} $outputfile]} {
        set outchan stdout
    } else {
        set outchan [open $outputfile "a"]
    }
    if {[string match {} $inputfile]} {
        puts $outchan "File: <stdin>"
    } else {
        puts $outchan "File: $inputfile"
    }
    if {![catch {file size $inputfile} filesize]} {
	puts $outchan "File size: $filesize bytes"
    }
    puts $outchan "Mesh title: $title"
    if {[string first "\n" $desc]<0} {
	# Single line description
	puts $outchan "Mesh  desc: $desc"
    } else {
	# Mult-line description
	regsub -all "\n" $desc "\n: " desc
	puts $outchan "Mesh desc---\n: $desc"
    }
    puts $outchan [GetMeshStructureInfo]
    if {![string match {} $outputfile]} {
        close $outchan
    }
    exit 0
}

if {[string match Vf_GeneralMesh3f [GetMeshType]] &&
	![string match irregular $gridmap($grid)]} {
    puts stderr "Input file has irregular grid; so must output."
}

if {!$quiet} {
    puts stderr "Copying Title field: $title\nYou may want to edit this\
	    field in the output by hand."
}

if {![string match "x:y:z" $flipstr] || ![string match {} $clipstr]} {
    # Perform coordinate transform and/or clipping
    CopyMesh 0 0 0 $flipstr $clipstr $cliprange
}

# $filename == "" --> write to stdout
WriteMesh $outputfile $formatmap($format) $gridmap($grid) $title $desc

exit 0
