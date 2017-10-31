
# @@ Meta Begin
# Package math 1.2.5
# Meta activestatetags ActiveTcl Public Tcllib
# Meta as::build::date 2013-08-13
# Meta as::origin      http://sourceforge.net/projects/tcllib
# Meta category        Tcl Math Library Tcl Math Library
# Meta description     Tcl Math Library Combinatorial functions in the Tcl
# Meta description     Math Library
# Meta license         BSD
# Meta platform        tcl
# Meta require         {Tcl 8.0}
# Meta subject         math statistics
# Meta summary         math math::combinatorics
# @@ Meta End


if {![package vsatisfies [package provide Tcl] 8.0]} return

package ifneeded math 1.2.5 [string map [list @ $dir] {
        # ACTIVESTATE TEAPOT-PKG BEGIN REQUIREMENTS

        package require Tcl 8.0

        # ACTIVESTATE TEAPOT-PKG END REQUIREMENTS

            source [file join {@} math.tcl]

        # ACTIVESTATE TEAPOT-PKG BEGIN DECLARE

        package provide math 1.2.5

        # ACTIVESTATE TEAPOT-PKG END DECLARE
    }]
