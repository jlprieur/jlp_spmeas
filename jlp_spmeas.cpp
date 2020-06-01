/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Name: jlp_spmeas.cpp SpmFrame class  
* Purpose: perform automatic measurements of speckle FITS autoc files
*          or theta calibration of FITS long integration files
*
* JLP 
* Version 13-03-2020
-------------------------------------------------------------------*/
#include <stdlib.h>   // exit()

/*
#define DEBUG
*/

#include "spm_frame.h"
// Definition of identifiers in "spm_frame_id.h"
#include "spm_frame_id.h"
#include "jlp_spmeas_utils.h" // JLP_Init_BIN_PARAM()

#if defined(__WXGTK__) || defined(__WXMOTIF__) || defined(__WXMAC__) || defined(__WXMGL__) || defined(__WXX11__)
    #define USE_XPM
#endif

#ifdef USE_XPM
    #include "mondrian.xpm"
#endif

BEGIN_EVENT_TABLE(SpmFrame, wxFrame)

// Menu Files:
   EVT_MENU(ID_LOAD_FITS_LIST, SpmFrame::OnLoadFITSList)
   EVT_MENU(ID_SAVE_TO_CSV,    SpmFrame::OnSaveToCsv)
   EVT_MENU(ID_SAVE_TO_LATEX,  SpmFrame::OnSaveToLatex)
   EVT_MENU(ID_SAVE_TO_TXT,    SpmFrame::OnSaveToTxt)
   EVT_MENU(ID_QUIT,           SpmFrame::OnQuit)

// Menu Config
   EVT_MENU(ID_SET_CALIBT,     SpmFrame::OnSetCalibTheta)
   EVT_MENU(ID_SET_AUTOM,      SpmFrame::OnSetAutoMeasures)
   EVT_MENU(ID_LOAD_AUTOMPARAM, SpmFrame::OnLoadAutomParamFile)

// Menu Start Run 
   EVT_MENU(ID_START_RUN,     SpmFrame::OnStartRun)

// Miscellaneous:
   EVT_MENU(ID_ABOUT,          SpmFrame::OnAbout)
   EVT_MENU(ID_HELP,           SpmFrame::OnHelp)

END_EVENT_TABLE()

//----------------------------------------------------------------------
// MyApp
//----------------------------------------------------------------------

class MyApp: public wxApp
{
public:
   bool OnInit();
};

IMPLEMENT_APP(MyApp)

bool MyApp::OnInit()
{
// Transform coma into point for numbers:
setlocale(LC_NUMERIC, "C");

    // use standard command line handling:
    if ( !wxApp::OnInit() )
        return false;

    // parse the cmd line
    int xx = 400, yy = 200;
    if ( argc == 3 )
    {
        wxSscanf(wxString(argv[1]), wxT("%d"), &xx);
        wxSscanf(wxString(argv[2]), wxT("%d"), &yy);
    }

// Create the main frame window
    SpmFrame *frame = new SpmFrame(_T("JLP_spmeas"), xx, yy);

// Give it an icon
// The wxICON() macros loads an icon from a resource under Windows
// and uses an #included XPM image under GTK+ and Motif
#ifdef USE_XPM
    frame->SetIcon( wxICON(mondrian) );
#endif

    frame->Show(true);

    return true;
}


/**********************************************************************
* SpmFrame constructor
*
* INPUT:
*   xx,yy : size of created window
*
***********************************************************************/
SpmFrame::SpmFrame(const wxChar *title, int xx, int yy)
       : wxFrame(NULL, wxID_ANY, title, wxPoint(-1, -1), wxSize(xx, yy))
{
wxString str1;
wxSize my_init_size;

m_StatusBar = NULL;

// Status bar:
// Create a status bar with two fields at the bottom:
  m_StatusBar = CreateStatusBar(2);
// First field has a variable length, second has a fixed length:
  int widths[2];
  widths[0] = -1;
  widths[1] = 200;
  SetStatusWidths( 2, widths );


// Min and Max size:
int iwidth, iheight;
  iwidth = xx;
  iheight = yy;
  my_init_size = wxSize(iwidth + 5, (iheight*6)/5 + 10);
  SetClientSize(my_init_size);

// Init private parameters:
  Spm_InitParameters();

// Create a menu on top of the window:
  Spm_SetupMenu();

initialized = 1234;

return;
}
/********************************************************************
* Init private parameters: 
********************************************************************/
void SpmFrame::Spm_InitParameters()
{
list_fname1[0] = '\0';
latex_fname1[0] = '\0';
param_fname1[0] = '\0';
csv_fname1[0] = '\0';
calib_fname1[0] = '\0';

// Processing mode: 1=auto-measurements 2=theta-calibration
pmode1 = 1;

// Load default values:
JLP_Init_BIN_PARAM(&bin_param1);

return;
}
/********************************************************************
* Setup the menu on top of main frame
********************************************************************/
void SpmFrame::Spm_SetupMenu()
{
SetHelpText( _T("Program to perform automatic measurements of speckle autoc files\n or theta calibration of long intergration files") );

  menu_bar = new wxMenuBar;

// ***************** File menu **********************************
  menuFile = new wxMenu;

  menuFile->Append(ID_LOAD_FITS_LIST, _T("Open a list of FITS files"),
                    _T("Open a list of FITS files"));
  menuFile->AppendSeparator();
  menuFile->Append( ID_SAVE_TO_CSV, _T("Select csv file for autom. measures"),
                    _T("Select csv file for automatic measures"));
  menuFile->Append( ID_SAVE_TO_LATEX, _T("Select latex file for autom. measures"),
                    _T("Select latex file for automatic measures"));
  menuFile->Append( ID_SAVE_TO_TXT, _T("Select output file for theta calib."),
                    _T("Select output file for theta calibration"));
  menuFile->AppendSeparator();

  menuFile->Append(ID_QUIT, _T("E&xit\tAlt-X"), _T("Quit program"));
  menu_bar->Append(menuFile, _T("&File"));

// ***************** Config menu **********************************
  menuConfig = new wxMenu;
  menuConfig->Append(ID_SET_CALIBT, _T("Set Theta Calibration Mode"),
                    _T("Set theta calibration mode"));
  menuConfig->Append(ID_SET_AUTOM, _T("Set Automatic Measurement Mode"),
                    _T("Set automatic measurement mode"));
  menuConfig->AppendSeparator();
  menuConfig->Append(ID_LOAD_AUTOMPARAM, _T("Load parameters for auto. measurements"),
                    _T("Load parameter file for automatic measurement mode"));
  menu_bar->Append(menuConfig, _T("&Configuration"));

// ***************** Start Run menu **********************************
  menuStartRun = new wxMenu;
  menuStartRun->Append(ID_START_RUN, _T("Start Run (with selected Mode)"),
                    _T("Start Run (with previously selected mode)"));
  menu_bar->Append(menuStartRun, _T("&Run"));

// ***************** Help menu ******************************
  menuHelp = new wxMenu;
  menuHelp->Append( ID_HELP, _T("Help"));
  menuHelp->Append( ID_ABOUT, _T("&About..."));
  menu_bar->Append(menuHelp, _T("Help"));

  SetMenuBar(menu_bar);

return;
}

void SpmFrame::OnQuit (wxCommandEvent& event)
{
// printf("OnQuit: calling Close\n");
    Close(true);
// event.Skip();
}

/*****************************************************************
* Help
*****************************************************************/
void SpmFrame::OnHelp( wxCommandEvent& WXUNUSED(event) )
{
 (void)wxMessageBox(_T("Sorry: \"Help\" is not implemented yet\n")
                    _T("Current version: May 2020"),
                    _T("JLP_SPMEAS"),
                     wxICON_INFORMATION | wxOK );
}
/*****************************************************************
* About
*****************************************************************/
void SpmFrame::OnAbout( wxCommandEvent& WXUNUSED(event) )
{
 (void)wxMessageBox( _T("JLP_SPMEAS\n")
                     _T("Jean-Louis Prieur (c) 2020\n")
                     _T("Created with wxWidgets"), _T("JLP_SPMEAS"),
                     wxICON_INFORMATION | wxOK );
}
/************************************************************************
* Display text in status bar
*************************************************************************/
void SpmFrame::SetText_to_StatusBar(wxString str1, const int icol)
{
// Update the first field (since 2nd argument is 0 here) of the status bar:
  if(m_StatusBar != NULL) m_StatusBar->SetStatusText(str1, icol);
}

