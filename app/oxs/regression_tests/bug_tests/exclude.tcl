lappend excludelist [list demagtest-c {3 4} "Bad numerics in demag kernel"]
lappend excludelist [list alphacheck 2 "Bad numerics?"]
lappend excludelist [list ovfout {3 4 5} "Bad numerics?"]
lappend excludelist [list fixedspins 7 "Bad numerics?"]
global tcl_platform
if {![info exists tcl_platform(wordSize)] || $tcl_platform(wordSize)<8} {
   lappend excludelist [list bigindex {} "Requires 64-bit OS and 128 GB of memory"]
} else {
   lappend excludelist [list bigindex {} "Requires 128 GB of memory (and many seconds)"]
}
