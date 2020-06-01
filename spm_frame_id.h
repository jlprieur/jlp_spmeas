/**************************************************************************
* Name: spm_frame_id.h  (SpmFrame class)
*
* Purpose: automatica measurements os speckle autoc files 
*
* list of ID used by SpmFrame class
*
* JLP
* Version 07/05/2020
****************************************************************************/
#ifndef _spm_frame_id_h_
#define _spm_frame_id_h_

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

//---------------------------------------------------------------------
enum
{
// Pb in menu found in 2015
//  ID_QUIT               = wxID_EXIT,
// Solution:
  ID_QUIT               = 999,
// File:
  ID_SAVE_TO_CSV = 1000,
  ID_SAVE_TO_LATEX,
  ID_SAVE_TO_TXT,
  ID_LOAD_FITS_LIST,

  ID_SET_CALIBT,
  ID_SET_AUTOM,
  ID_LOAD_AUTOMPARAM,
  ID_START_RUN,

// Help:
  ID_ABOUT             = wxID_ABOUT,
  ID_HELP              = wxID_HELP
};

#endif

