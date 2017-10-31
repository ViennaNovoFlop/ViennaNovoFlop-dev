/* FILE: chanwrap.cc                 -*-Mode: c++-*-
 *
 * C++ wrappers for Tcl channels.  These wrappers are exception safe,
 * and also provide low-overhead, buffered byte-at-a-time access to
 * the channels.  Support is only for strictly read-only or write-only
 * channels.  Also, once a channel is attached to one of these objects,
 * only access through the object is allowed; e.g., direct calls to
 * Tcl_Read will almost certainly result in undetected reordering or
 * lossage of bytes on the input channel.  If the channel supports
 * random seeks (via Tcl_Tell and Tcl_Seek), then the channel can
 * be detached from the object, after which all the standard Tcl
 * I/O calls can be used again.
 *
 * NOTICE: Please see the file ../../LICENSE
 *
 * Last modified on: $Date: 2007-03-21 23:29:34 $
 * Last modified by: $Author: donahue $
 */

#include <assert.h>

#include "oc.h"

#include "chanwrap.h"
#include "functions.h"

/* End includes */     /* Optional directive to build.tcl */


////////////////////////////////////////////////////////////////////////
// Wrapper for read-access Tcl channels

#define NB_INPUTCHANNEL_BUFSIZE 65536
#define NB_INPUTCHANNEL_ERRBUFSIZE 4096

Nb_InputChannel::Nb_InputChannel()
  : chan(NULL),buf_fillsize(0),buf_offset(0)
{
    assert(NB_INPUTCHANNEL_BUFSIZE>0);
    buf.SetLength(NB_INPUTCHANNEL_BUFSIZE);
}

void Nb_InputChannel::OpenFile(const char* filename,
                               const char* translation)
{
  Close();

  Oc_DirectPathname(filename,file_descriptor);

  chan = Tcl_OpenFileChannel(NULL,file_descriptor.GetStr(),
			     OC_CONST84_CHAR("r"),0);
  if(chan==NULL) {
    // Unable to open file.
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"Nb_InputChannel","OpenFile",
                       NB_INPUTCHANNEL_ERRBUFSIZE+2000,
                       "Unable to open file \"%.2000s\" for reading.",
                       file_descriptor.GetStr()));
  }

  if(Tcl_SetChannelOption(NULL,chan,
			  OC_CONST84_CHAR("-translation"),
			  OC_CONST84_CHAR(translation))!=TCL_OK ||
     Tcl_SetChannelOption(NULL,chan,
			  OC_CONST84_CHAR("-buffering"),
			  OC_CONST84_CHAR("full"))!=TCL_OK ||
     Tcl_SetChannelOption(NULL,chan,OC_CONST84_CHAR("-buffersize"),
         OC_CONST84_CHAR(OC_STRINGIFY(NB_INPUTCHANNEL_BUFSIZE)))!=TCL_OK) {
    // Channel configuration error.
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"Nb_InputChannel","OpenFile",
                       NB_INPUTCHANNEL_ERRBUFSIZE+2000,
                       "Unable to configure input channel"
                       " for file \"%.2000s\".",
                       file_descriptor.GetStr()));
  }
  buf_fillsize = buf_offset = 0;
}

int Nb_InputChannel::FillBuf()
{
  buf_offset=0;
  buf_fillsize = Tcl_Read(chan,(char*)buf,NB_INPUTCHANNEL_BUFSIZE);
  if(buf_fillsize<0) {
    buf_fillsize = 0;
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"Nb_InputChannel","FillBuf",
                       NB_INPUTCHANNEL_ERRBUFSIZE+4000,
                       "Read error on \"%.2000s\": %.2000s",
                       file_descriptor.GetStr(),
                       Tcl_ErrnoMsg(Tcl_GetErrno())));
  }
  return buf_fillsize;
}

#undef NB_INPUTCHANNEL_BUFSIZE
#undef NB_INPUTCHANNEL_ERRBUFSIZE

////////////////////////////////////////////////////////////////////////
// Wrapper for write-access Tcl channels

#define NB_OUTPUTCHANNEL_ERRBUFSIZE 4096

Nb_OutputChannel::Nb_OutputChannel()
  : chan(NULL),buf_offset(0)
{
  assert(NB_OUTPUTCHANNEL_BUFSIZE>0);
  buf.SetLength(NB_OUTPUTCHANNEL_BUFSIZE);
}

void Nb_OutputChannel::OpenFile(const char* filename,
                               const char* translation)
{
  Close();

  Oc_DirectPathname(filename,file_descriptor);

  chan = Tcl_OpenFileChannel(NULL,file_descriptor.GetStr(),
			     OC_CONST84_CHAR("w"),0);
  if(chan==NULL) {
    // Unable to open file.
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"Nb_OutputChannel","OpenFile",
                       NB_OUTPUTCHANNEL_ERRBUFSIZE+2000,
                       "Unable to open file \"%.2000s\" for writing.",
                       file_descriptor.GetStr()));
  }

  if(Tcl_SetChannelOption(NULL,chan,OC_CONST84_CHAR("-translation"),
			  OC_CONST84_CHAR(translation))!=TCL_OK ||
     Tcl_SetChannelOption(NULL,chan,OC_CONST84_CHAR("-buffering"),
			  OC_CONST84_CHAR("full"))!=TCL_OK ||
     Tcl_SetChannelOption(NULL,chan,OC_CONST84_CHAR("-buffersize"),
         OC_CONST84_CHAR(OC_STRINGIFY(NB_OUTPUTCHANNEL_BUFSIZE)))!=TCL_OK) {
    // Channel configuration error.
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"Nb_OutputChannel","OpenFile",
                       NB_OUTPUTCHANNEL_ERRBUFSIZE+2000,
                       "Unable to configure output channel"
                       " for file \"%.2000s\".",
                       file_descriptor.GetStr()));
  }
  buf_offset = 0;
}

int Nb_OutputChannel::WriteBuf()
{
  int bytes_written = Tcl_Write(chan,(char*)buf,buf_offset);
  if(bytes_written<0) {
    OC_THROW(Oc_Exception(__FILE__,__LINE__,"Nb_OutputChannel","WriteBuf",
                       NB_OUTPUTCHANNEL_ERRBUFSIZE+4000,
                       "Write error on \"%.2000s\": %.2000s",
                       file_descriptor.GetStr(),
                       Tcl_ErrnoMsg(Tcl_GetErrno())));
  }
  buf_offset=0;
  return bytes_written;
}

#undef NB_OUTPUTCHANNEL_BUFSIZE
#undef NB_OUTPUTCHANNEL_ERRBUFSIZE
