if {[catch {package present Tcl 8.5.0}]} { return }
package ifneeded Tk 8.5.12 [list load [file join $dir .. libtk8.5.so] Tk]
