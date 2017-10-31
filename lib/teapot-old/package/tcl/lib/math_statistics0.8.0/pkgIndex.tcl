
# @@ Meta Begin
# Package math::statistics 0.8.0
# Meta activestatetags ActiveTcl Public Tcllib
# Meta as::build::date 2013-08-13
# Meta as::origin      http://sourceforge.net/projects/tcllib
# Meta category        Tcl Math Library
# Meta description     Basic statistical functions and procedures
# Meta license         BSD
# Meta platform        tcl
# Meta require         {Tcl 8.4}
# Meta require         math
# Meta require         {math::linearalgebra 1.0}
# Meta subject         mathematics {data analysis} statistics
# Meta summary         math::statistics
# @@ Meta End


if {![package vsatisfies [package provide Tcl] 8.4]} return

package ifneeded math::statistics 0.8.0 [string map [list @ $dir] {
        # ACTIVESTATE TEAPOT-PKG BEGIN REQUIREMENTS

        package require Tcl 8.4
        package require math
        package require math::linearalgebra 1.0

        # ACTIVESTATE TEAPOT-PKG END REQUIREMENTS

            source [file join {@} statistics.tcl]

        # ACTIVESTATE TEAPOT-PKG BEGIN DECLARE

        package provide math::statistics 0.8.0

        # ACTIVESTATE TEAPOT-PKG END DECLARE
    }]
