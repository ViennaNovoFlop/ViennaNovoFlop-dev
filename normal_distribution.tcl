package require math::statistics
set n 100000
set m  1.
set s  1.

set l [list 1 2 3]

proc mittel { dat } {
set i 0
set k 0.
foreach a $dat {
incr i
set k [expr {$k+$a}]
#puts $a
}
set res [expr {$k/$i}]
return $res
}

#set l [::math::satistics::random-normal $m $s $n]
#puts  [mittel [::math::statistics::random-normal $m $s $n]]
#puts  [::math::statistics::mean [::math::statistics::random-normal $m $s $n]]
#puts $l
#puts [mittel $l]
# now create list for x, y and z component each direction normal distributed
# ten times three vectors
for {set x 0} {$x < $n} {incr x} {
lappend vec [::math::statistics::random-normal 0. 1. 1] 
}
#puts $vec
puts "mean:  [::math::statistics::mean $vec]"
puts "stdev: [::math::statistics::stdev $vec]"

