/* File: crc.h
 *
 * Table-based 32-bit crc calculator.
 *
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2010-07-16 22:33:54 $
 * Last modified by: $Author: donahue $
 */

#ifndef _NB_CRC
#define _NB_CRC

#include "oc.h"

/* End includes */     /* Optional directive to pimake */

OC_UINT4m Nb_ComputeCRC(OC_UINT4m size,const unsigned char* bytebuf);
OC_UINT4m Nb_ComputeCRC(Tcl_Channel chan,OC_UINT4m* bytes_read=NULL);

// Tcl wrappers for the above.  The first returns the CRC,
// the second returns a two element list consisting of the CRC
// and the number of bytes read.
Tcl_CmdProc NbComputeCRCBufferCmd;
Tcl_CmdProc NbComputeCRCChannelCmd;

#endif // _NB_CRC
