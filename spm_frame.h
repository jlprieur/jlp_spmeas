/********************************************************
* Name: spm_frame.h SpmFrame class
* Purpose: perform automatic measurements of speckle FITS autoc files
*          or theta calibration of FITS long integration files
*
* JLP
* Version 13/05/2020 
*********************************************************/
#ifndef _spm_frame__
#define _spm_frame__

#include <stdlib.h> // exit()

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "jlp_spmeas_def.h"

//----------------------------------------------------------------------
// class definitions
//----------------------------------------------------------------------

class SpmFrame: public wxFrame
{
public:

// In "jlp_spmeas.cpp":
    SpmFrame(const wxChar *title, int x, int y);
// Not exiting...
//    ~SpmFrame() {printf("DDestructor called\n"); Destroy(); return;};
    ~SpmFrame() {
     // printf("exit function called\n"); 
      exit(-1);
      };
   
    void Spm_InitParameters();
    void Spm_SetupMenu();
    void OnQuit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
    void OnHelp(wxCommandEvent& event);

    void SetText_to_StatusBar(wxString text, const int icol);
 
// In "spm_frame_menu.cpp":
    void OnSaveToCsv(wxCommandEvent& event);
    void OnSaveToLatex(wxCommandEvent& event);
    void OnSaveToTxt(wxCommandEvent& event);
    void OnLoadFITSList(wxCommandEvent &event);
    int  LoadFITSList();

    void OnSetCalibTheta(wxCommandEvent& event);
    void OnSetAutoMeasures(wxCommandEvent& event);
    void OnLoadAutomParamFile(wxCommandEvent& event);
    int LoadAutomParamFile();
    void OnStartRun(wxCommandEvent& event);
    int Update_BIN_PARAM();

// In "spm_frame_autom.cpp":
    int ProcessListOfAutocFilesForAutom();
    int AutoMeasureOfAutoc(char *in_name);

// In "spm_frame_theta_calib.cpp":
    int ProcessListForThetaCalib();

private:
  wxStatusBar *m_StatusBar;
  int nx1, ny1, initialized, pmode1;
  BIN_PARAM bin_param1;

  char list_fname1[256], latex_fname1[256];
  char csv_fname1[256], param_fname1[256], calib_fname1[256];

// Menus:
  wxMenuBar *menu_bar;
  wxMenu *menuFile, *menuHelp;
  wxMenu *menuConfig, *menuStartRun;

  DECLARE_EVENT_TABLE()

};

#endif
