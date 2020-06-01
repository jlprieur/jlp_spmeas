/******************************************************************************
* Name:        spm_frame_menu (SpmFrame class)
* Purpose:     handling menu events of SpmFrame class (see jlp_spmeas.cpp)
* Author:      JLP
* Version:     07/05/2020
******************************************************************************/
#include "spm_frame.h"
#include "jlp_spmeas_dlg.h"
#include "jlp_spmeas_utils.h"

/* Declared in spm_frame.h

void OnSaveToCsv(wxCommandEvent &WXUNUSED(event));
void OnLoadFITSList(wxCommandEvent &WXUNUSED(event));
int LoadFITSList();

*/

/************************************************************************
* Select a csv (comma separated value) file for output
************************************************************************/
void SpmFrame::OnSaveToCsv(wxCommandEvent &WXUNUSED(event))
{
wxString filename, str1;
FILE *fp_out;

if(initialized != 1234) return;

 csv_fname1[0] = '\0';
 wxFileDialog dialog(NULL, _T("Select cvx text file for output"), wxEmptyString,
                     wxEmptyString, _T("Files (*.txt)|*.txt"),
                     wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
 if (dialog.ShowModal() != wxID_OK) return;

 filename = dialog.GetPath();
 strcpy(csv_fname1, filename.mb_str());

// Try to open it:
if((fp_out = fopen(csv_fname1,"w")) == NULL) {
  str1.Printf(_T("Error opening output file >%s<\n"), csv_fname1);
  wxMessageBox(str1, _T("OnSaveToCsv"), wxICON_ERROR);
  csv_fname1[0] = '\0';
  return;
  }

fclose(fp_out);
return;
}
/************************************************************************
* Select a txt file for output
************************************************************************/
void SpmFrame::OnSaveToTxt(wxCommandEvent &WXUNUSED(event))
{
wxString filename, str1;
FILE *fp_out;

if(initialized != 1234) return;

 calib_fname1[0] = '\0';
 wxFileDialog dialog(NULL, _T("Select txt file for calib. output"), 
                     wxEmptyString, wxEmptyString, _T("Files (*.txt)|*.txt"),
                     wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
 if (dialog.ShowModal() != wxID_OK) return;

 filename = dialog.GetPath();
 strcpy(calib_fname1, filename.mb_str());

// Try to open it:
if((fp_out = fopen(calib_fname1,"w")) == NULL) {
  str1.Printf(_T("Error opening output file >%s<\n"), calib_fname1);
  wxMessageBox(str1, _T("OnSaveToTxt"), wxICON_ERROR);
  calib_fname1[0] = '\0';
  return;
  }

fclose(fp_out);
return;
}
/************************************************************************
* Select a latex file for output
************************************************************************/
void SpmFrame::OnSaveToLatex(wxCommandEvent &WXUNUSED(event))
{
wxString filename, str1;
FILE *fp_out;

if(initialized != 1234) return;

 latex_fname1[0] = '\0';
 wxFileDialog dialog(NULL, _T("Select latex file for output"), wxEmptyString,
                     wxEmptyString, _T("Files (*.tex)|*.tex"),
                     wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
 if (dialog.ShowModal() != wxID_OK) return;

 filename = dialog.GetPath();
 strcpy(latex_fname1, filename.mb_str());

// Try to open it:
if((fp_out = fopen(latex_fname1,"w")) == NULL) {
  str1.Printf(_T("Error opening output file >%s<\n"), latex_fname1);
  wxMessageBox(str1, _T("OnSaveToLatex"), wxICON_ERROR);
  csv_fname1[0] = '\0';
  return;
  }

fclose(fp_out);
return;
}
/******************************************************************
* Set Calib Theta mode
*******************************************************************/
void SpmFrame::OnSetCalibTheta(wxCommandEvent& event)
{
// Processing mode: 1=auto-measurements 2=theta-calibration
pmode1 = 2;
return;
}
/******************************************************************
* Set Auto Measures mode
*******************************************************************/
void SpmFrame::OnSetAutoMeasures(wxCommandEvent& event)
{
// Processing mode: 1=auto-measurements 2=theta-calibration
pmode1 = 1;
return;
}
/******************************************************************
* Load parameter file for the Auto Measures mode
*******************************************************************/
void SpmFrame::OnLoadAutomParamFile(wxCommandEvent& event)
{
if(initialized != 1234) return;

// Processing mode: 1=auto-measurements 2=theta-calibration
if(pmode1 != 1) {
 wxMessageBox(_T("Error: please select first auto-measurement mode!"), 
                 _T("OnLoadAutomParamFile"), wxICON_ERROR);
 return;
 }

LoadAutomParamFile();

return;
}
/******************************************************************
* Load parameter file for the Auto Measures mode
*******************************************************************/
int SpmFrame::LoadAutomParamFile()
{
wxString param_filename, str1;
int status;
FILE *fp_in;

if(initialized != 1234) return(-1);

param_fname1[0] = '\0';

// Prompt the user for the parameter filename
param_filename = wxFileSelector(_T("Select parameter file for automatic measurements"),
                 _T(""), _T(""), _T("txt"), _T("txt files (*.txt)|*.txt"));

if ( param_filename.empty() ) return(-1);

strcpy(param_fname1, (const char*)param_filename.mb_str());

// Read it and load it to private variable:
JLP_Load_BIN_PARAM_from_file(param_fname1, &bin_param1);

return(0);
}
/******************************************************************
* Start run 
*******************************************************************/
void SpmFrame::OnStartRun(wxCommandEvent& event)
{
int status;

  if(initialized != 1234) return;

  switch (pmode1) { 
     default:
     case 1:
       status = Update_BIN_PARAM();
       if(status == 0) ProcessListOfAutocFilesForAutom();
       break;
     case 2:
       ProcessListForThetaCalib();
       break;
     }
return;
}
/******************************************************************
* Start run
*******************************************************************/
int SpmFrame::Update_BIN_PARAM()
{
int dlg_answer, status = -1;
JLP_Spmeas_Dlg *SpmeasDlg;

// Load parameters from file if necessary:
if(param_fname1[0] != '\0')
   JLP_Load_BIN_PARAM_from_file(param_fname1, &bin_param1);

// Check whether the user agrees:
SpmeasDlg = new JLP_Spmeas_Dlg(NULL, bin_param1, 
                               wxT("Parameters for automatic measurements"));

// WARNING: There is a bug here coming from "wxwidgets"
// when using ShowModal, the system doesn't return
// The computer may even crash if there are too many of those processed hanging
// around and using CPU time !
 dlg_answer = SpmeasDlg->ShowModal();

// If OK, exit from loop:
 if(dlg_answer == 0) {
   SpmeasDlg->RetrieveData(&bin_param1);
   status = 0;
 }

delete SpmeasDlg;

return(status);
}
/******************************************************************
* Load a new FITS list file
*******************************************************************/
void SpmFrame::OnLoadFITSList(wxCommandEvent &WXUNUSED(event))
{

  if(initialized != 1234) return;

  LoadFITSList();

return;
}
/******************************************************************
* Load a list of 2D FITS files
*
*******************************************************************/
int SpmFrame::LoadFITSList()
{
wxString list_filename, str1;
int status;
FILE *fp_in;

if(initialized != 1234) return(-1);

  list_fname1[0] = '\0';
// Prompt the user for the FITS list filename
  list_filename = wxFileSelector(_T("Select image FITS list file"), 
                        _T(""), _T(""), _T("txt"),
                      _T("txt files (*.txt)|*.txt"));

  if ( list_filename.empty() ) return(-1);

strcpy(list_fname1, (const char*)list_filename.mb_str());

// Try to open it:
if((fp_in = fopen(list_fname1,"r")) == NULL) {
  str1.Printf(_T("Error opening input file >%s<\n"), list_fname1);
  wxMessageBox(str1, _T("LoadFITSList"), wxICON_ERROR);
  list_fname1[0] = '\0';
  return(-1);
  }

fclose(fp_in);
return(0);
}
