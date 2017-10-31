#!/bin/sh
#
# Bourne shell script to convert *.h and *.cc files using
#
#           INT*/UINT*/REAL*/BOOL/...
#
# typedefs to
#
#        OC_INT*/OC_UINT*/OC_REAL*/OC_BOOL/...
#
# typedefs.  This change was made in the main OOMMF sources
# on 16-Jul-2010.  This script may be useful to OOMMF extension
# authors to convert sources from OOMMF 1.2a3 to OOMMF 1.2a4.

for i in `find . -xdev -type f \( -name \*.h -o -name \*.cc \) -print` ; do echo "---- $i ---" ; cat $i | sed 's/\(^\|[^_U]\)\(\(U\|\)INT\(2\|4\|8\)\)/\1OC_\2/g' > foo ; mv -f foo $i ; done

 for i in `find . -xdev -type f \( -name \*.h -o -name \*.cc \) -print` ; do echo "---- $i ---" ; cat $i | sed 's/\(^\|[^_WTUS]\)\(BOOL\|BYTE\|\(U\|S\|\)CHAR\)/\1OC_\2/g' > foo ; mv -f foo $i ; done

for i in `find . -xdev -type f \( -name \*.h -o -name \*.cc \) -print` ; do echo "---- $i ---" ; cat $i | sed 's/\(^\|[^_]\)\(\(SQRT_\|CUBE_ROOT_\|\)REAL\(4\|8\|WIDE\)\)/\1OC_\2/g' > foo ; mv -f foo $i ; done

for i in `find . -xdev -type f \( -name \*.h -o -name \*.cc \) -print` ; do echo "---- $i ---" ; cat $i | sed 's/\(^\|[^_]\)\(\(POINTERWIDTH\|FP_REGISTER_EXTRA_PRECISION\)\)/\1OC_\2/g' > foo ; mv -f foo $i ; done

