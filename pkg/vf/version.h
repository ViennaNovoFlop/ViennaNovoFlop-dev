/* FILE: version.h                    -*-Mode: c++-*-
 *
 *	Version macros for the OOMMF Vector Field extension.
 *
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2004-09-01 15:48:49 $
 * Last modified by: $Author: dgp $
 */

#ifndef _VF_VERSION
#define _VF_VERSION

#include "oc.h"

/* End includes */

/*
 * NOTE: When version number information is changed here, it must
 * also be changed in the following files:
 *	./vf.tcl
 */

#define VF_MAJOR_VERSION	1
#define VF_MINOR_VERSION	2
#define VF_RELEASE_LEVEL	0
#define VF_RELEASE_SERIAL	4

#define VF_VERSION OC_MAKE_VERSION(VF)

#endif /* _VF_VERSION */
