#!/bin/sh
# Shell wrapper around the main installation script.  Will use the
# wish or tclsh inside of the distribution to execute this
# installer.
#
# Auto-detects the location of the distribution. Auto-detects
# availability of X and chooses between gui and terminal based
# installers using this information.

dist=`dirname $0`

# initialize the environment so that the distributed "wish" and
# "tclsh" are able to find their libraries despite being crippled with
# the special path value we will replace later during the installation
# with the actual path to the installation.

TCL_LIBRARY=$dist/payload/lib/tcl8.5
TK_LIBRARY=$dist/payload/lib/tk8.5
LD_LIBRARY_PATH=$dist/payload/lib
DYLD_LIBRARY_PATH=$dist/payload/lib
SHLIB_PATH=$dist/payload/lib
LIBPATH=$dist/payload/lib

export TCL_LIBRARY TK_LIBRARY LD_LIBRARY_PATH DYLD_LIBRARY_PATH SHLIB_PATH LIBPATH

# Determine availability of X and choose an installer based on that
# information.

if tty -s; then                 # Do we have a terminal?
    if [ x"$DISPLAY" != x -a x"$xterm_loop" = x ]; then  # No, but do we have X?

	# Check for valid DISPLAY variable

	if [ `echo exit | $dist/payload/bin/wish8.5 2>&1 | grep fail | wc -l` -eq 0 ]
	then
	    echo _____________________________________________
	    echo Launching graphical installer on $DISPLAY
	    echo ...

	    $dist/payload/bin/wish8.5 $dist/install.tk "$@"

	    # pwd = inside the unpacked distribution ...
	    # go one level up and remove the directory
	    #cd ..
	    #rm -rf $dist
	    exit
	fi
    fi
fi

# No X, use the terminal based installer 

$dist/payload/bin/tclsh8.5 $dist/install.tcl "$@"

# pwd = inside the unpacked distribution ...
# go one level up and remove the directory
#cd ..
#rm -rf $dist
exit
