/******************************************************************************
* jlp_spmeas_dlg.cpp
* Dialog box used for selecting the binary parameters for automatic meaurements 
*
* Author:      JLP 
* Version:     18/05/2020
******************************************************************************/
#include "jlp_spmeas_dlg.h"
#include "jlp_spmeas_utils.h" // JLP_Init_BIN_PARAM()

//*************************************************************************
enum
{
   ID_BPARAM_RELERRMAX = 800,
   ID_BPARAM_ERTMAX,
   ID_BPARAM_RMIN,
   ID_BPARAM_NLRMIN,
   ID_BPARAM_NLRMAX,
   ID_BPARAM_NEGMAX,
   ID_BPARAM_OK,
   ID_BPARAM_CANCEL,
   ID_BPARAM_DEFAULT
};

BEGIN_EVENT_TABLE(JLP_Spmeas_Dlg, wxDialog)
EVT_BUTTON  (ID_BPARAM_OK, JLP_Spmeas_Dlg::OnOKButton)
EVT_BUTTON  (ID_BPARAM_CANCEL, JLP_Spmeas_Dlg::OnCancelButton)
EVT_BUTTON  (ID_BPARAM_DEFAULT, JLP_Spmeas_Dlg::OnDefaultButton)
EVT_TEXT    (ID_BPARAM_RELERRMAX, JLP_Spmeas_Dlg::OnChangeParam)
EVT_TEXT    (ID_BPARAM_ERTMAX, JLP_Spmeas_Dlg::OnChangeParam)
EVT_TEXT    (ID_BPARAM_RMIN, JLP_Spmeas_Dlg::OnChangeParam)
EVT_TEXT    (ID_BPARAM_NLRMIN, JLP_Spmeas_Dlg::OnChangeParam)
EVT_TEXT    (ID_BPARAM_NLRMAX, JLP_Spmeas_Dlg::OnChangeParam)
EVT_TEXT    (ID_BPARAM_NEGMAX, JLP_Spmeas_Dlg::OnChangeParam)
END_EVENT_TABLE()

/********************************************************************
* Constructor:
*
********************************************************************/
JLP_Spmeas_Dlg::JLP_Spmeas_Dlg(wxFrame *parent, BIN_PARAM bin_param0,
                               const wxString title) 
        : wxDialog(parent, -1, title, wxPoint(400,100), wxDefaultSize,
                   wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER)
{

// Save input parameters to private variables:
bin_param1.error_theta10_deg_max = bin_param0.error_theta10_deg_max;
bin_param1.rel_error_rho10_max = bin_param0.rel_error_rho10_max;
bin_param1.rho10_min = bin_param0.rho10_min;
bin_param1.n_left_over_right_min = bin_param0.n_left_over_right_min;
bin_param1.n_left_over_right_max = bin_param0.n_left_over_right_max;
bin_param1.neg_max_for_patch = bin_param0.neg_max_for_patch;

initialized = 0;

// To avoid initialization problems with Windows:
// (An event is sent to "ChangeParameters"
//  as soon as a Text control is created...)
// First initialize the pointers to NULL
  TextCtrl_ertmax = NULL;
  TextCtrl_relerrmax = NULL;
  TextCtrl_rmin = NULL;
  TextCtrl_nlrmin = NULL;
  TextCtrl_nlrmax = NULL;
  TextCtrl_negmax = NULL;

  InitPanel();

initialized = 1234;

return;
}
/********************************************************************
* Init main panel 
*
********************************************************************/
void JLP_Spmeas_Dlg::InitPanel()
{
wxBoxSizer *topsizer, *button_sizer;
wxStaticBoxSizer *param_sizer;
int i, irows, icols, vgap = 12, hgap = 12, wwidth = 160;
wxFlexGridSizer *fgs1;
wxString ertmax_str, relerrmax_str;
wxString rmin_str, nlrmin_str, nlrmax_str, negmax_str;

// Initialize the text strings:
  ertmax_str.Printf(_T("%.2f"), bin_param1.error_theta10_deg_max);
  relerrmax_str.Printf(_T("%.2f"), bin_param1.rel_error_rho10_max);
  rmin_str.Printf(_T("%.2f"), bin_param1.rho10_min);
  nlrmin_str.Printf(_T("%.2f"), bin_param1.n_left_over_right_min);
  nlrmax_str.Printf(_T("%.2f"), bin_param1.n_left_over_right_max);
  negmax_str.Printf(_T("%.2f"), bin_param1.neg_max_for_patch);

// Create the text controls: 
  TextCtrl_relerrmax = new wxTextCtrl(this, ID_BPARAM_RELERRMAX, relerrmax_str);
  TextCtrl_ertmax = new wxTextCtrl(this, ID_BPARAM_ERTMAX, ertmax_str);
  TextCtrl_rmin = new wxTextCtrl(this, ID_BPARAM_RMIN, rmin_str);
  TextCtrl_nlrmin = new wxTextCtrl(this, ID_BPARAM_NLRMIN, nlrmin_str);
  TextCtrl_nlrmax = new wxTextCtrl(this, ID_BPARAM_NLRMAX, nlrmax_str);
  TextCtrl_negmax = new wxTextCtrl(this, ID_BPARAM_NEGMAX, negmax_str);

topsizer = new wxBoxSizer( wxVERTICAL );

/**********************************************************/
// Sizer surrounded with a rectangle, with a title on top:
 param_sizer = new wxStaticBoxSizer(wxVERTICAL, this, _T("Select icheck parameters"));
 irows = 6;
 icols = 2;
 wwidth = 100;
 fgs1 = new wxFlexGridSizer(irows, icols, vgap, hgap);
 fgs1->Add( new wxStaticText( this, wxID_ANY, 
                                    _T("1. relative error rho10 max :") ),
                  0, wxALIGN_LEFT|wxALIGN_CENTRE_VERTICAL|wxRIGHT, 10); 
 fgs1->Add(TextCtrl_relerrmax, 0, wxRIGHT, 10); 
 fgs1->Add( new wxStaticText( this, wxID_ANY, 
                                    _T("2. error theta10 max (deg):") ),
                  0, wxALIGN_LEFT|wxALIGN_CENTRE_VERTICAL|wxRIGHT, 10); 
 fgs1->Add(TextCtrl_ertmax, 0, wxRIGHT, 10); 
 fgs1->Add( new wxStaticText( this, wxID_ANY, 
                                    _T("3. rho10 min (pixels):") ),
                  0, wxALIGN_LEFT|wxALIGN_CENTRE_VERTICAL|wxRIGHT, 10); 
 fgs1->Add(TextCtrl_rmin, 0, wxRIGHT, 10); 
 fgs1->Add( new wxStaticText( this, wxID_ANY, 
                                    _T("4. n left over right min:") ),
                  0, wxALIGN_LEFT|wxALIGN_CENTRE_VERTICAL|wxRIGHT, 10); 
 fgs1->Add(TextCtrl_nlrmin, 0, wxRIGHT, 10); 
 fgs1->Add( new wxStaticText( this, wxID_ANY, 
                                    _T("5. n left over right max:") ),
                  0, wxALIGN_LEFT|wxALIGN_CENTRE_VERTICAL|wxRIGHT, 10); 
 fgs1->Add(TextCtrl_nlrmax, 0, wxRIGHT, 10); 
 fgs1->Add( new wxStaticText( this, wxID_ANY, 
                                    _T("6. neg. percent max for patch:") ),
                  0, wxALIGN_LEFT|wxALIGN_CENTRE_VERTICAL|wxRIGHT, 10); 
 fgs1->Add(TextCtrl_negmax, 0, wxRIGHT, 10); 

 param_sizer->Add(fgs1, 0, wxALL, 20);
 topsizer->Add(param_sizer, 0, wxALIGN_CENTER);

/**********************************************************/
button_sizer = new wxBoxSizer( wxHORIZONTAL );

//create two buttons that are horizontally unstretchable, 
  // with an all-around border with a width of 10 and implicit top alignment
 button_sizer->Add(
    new wxButton(this, ID_BPARAM_OK, _T("OK") ), 0, wxALL, 10);

 button_sizer->Add(
   new wxButton(this, ID_BPARAM_CANCEL, _T("Cancel") ), 0, wxALL, 10);

 button_sizer->Add(
   new wxButton(this, ID_BPARAM_DEFAULT, _T("Default values") ), 0, wxALL, 10);

  //create a sizer with no border and centered horizontally
 topsizer->Add(button_sizer, 0, wxALIGN_CENTER);

 SetSizer(topsizer);      // use the sizer for layout

 topsizer->SetSizeHints( this );   // set size hints to honour minimum size

return;
}
/**************************************************************************
* Handle "Cancel" button:
**************************************************************************/
void JLP_Spmeas_Dlg::OnCancelButton( wxCommandEvent& WXUNUSED(event) )
{
 if(initialized != 1234) return;

// Close dialog and return status = 1:
  EndModal(1); 
}
/**************************************************************************
* Handle "OK" button:
**************************************************************************/
void JLP_Spmeas_Dlg::OnOKButton( wxCommandEvent& WXUNUSED(event) )
{
double value;

 if(initialized != 1234) return;

// First check that fields contain real/integer values
// (and hence that have been loaded to the private variables...)
    if(!TextCtrl_ertmax->GetValue().ToDouble(&value)
       || !TextCtrl_relerrmax->GetValue().ToDouble(&value)
       || !TextCtrl_rmin->GetValue().ToDouble(&value)
       || !TextCtrl_nlrmin->GetValue().ToDouble(&value)
       || !TextCtrl_nlrmax->GetValue().ToDouble(&value)) {
        wxMessageBox(_T("Bad values: please first correct the fields!"),
                     _T("Spmeas parameters"), wxICON_ERROR);
        return;
       }

  if(DataIsOK() != true) {
// If bad new value restore previous value:
     wxMessageBox(_T("Bad parameters! (check that values are coherent or cancel)"),
                  _T(""), wxICON_ERROR);
   } else {
// Close dialog and return status = 0:
     EndModal(0);
   }

}
/**************************************************************************
* Handle "Default values" button:
**************************************************************************/
void JLP_Spmeas_Dlg::OnDefaultButton( wxCommandEvent& WXUNUSED(event) )
{
double value;

 if(initialized != 1234) return;

wxString ertmax_str, relerrmax_str;
wxString rmin_str, nlrmin_str, nlrmax_str, negmax_str;

// Load Default values:
  JLP_Init_BIN_PARAM(&bin_param1); 

// Initialize the text strings:
  ertmax_str.Printf(_T("%.2f"), bin_param1.error_theta10_deg_max);
  relerrmax_str.Printf(_T("%.2f"), bin_param1.rel_error_rho10_max);
  rmin_str.Printf(_T("%.2f"), bin_param1.rho10_min);
  nlrmin_str.Printf(_T("%.2f"), bin_param1.n_left_over_right_min);
  nlrmax_str.Printf(_T("%.2f"), bin_param1.n_left_over_right_max);
  negmax_str.Printf(_T("%.2f"), bin_param1.neg_max_for_patch);

   TextCtrl_ertmax->SetValue(ertmax_str);
   TextCtrl_relerrmax->SetValue(relerrmax_str);
   TextCtrl_rmin->SetValue(rmin_str);
   TextCtrl_nlrmin->SetValue(nlrmin_str);
   TextCtrl_nlrmax->SetValue(nlrmax_str);
   TextCtrl_negmax->SetValue(negmax_str);

// Refresh screen:
   Refresh();
return;
}
/**************************************************************************
* Handle text editing 
*
    TextCtrl_x
    TextCtrl_y
    TextCtrl_radius
    TextCtrl_radius_fact
    TextCtrl_poly_order
    TextCtrl_noise
*
**************************************************************************/
void JLP_Spmeas_Dlg::OnChangeParam( wxCommandEvent& event )
{
double old_value, new_value;
int old_ivalue;
long new_ivalue;
int status = 0;
wxString w_str;

if(initialized != 1234) return;

// First check that all text controls are created:
 if(!TextCtrl_ertmax || !TextCtrl_relerrmax
    || !TextCtrl_rmin || !TextCtrl_nlrmin || !TextCtrl_nlrmax) return; 

  switch (event.GetId())
  {
   case ID_BPARAM_RELERRMAX:
    {
    if(TextCtrl_relerrmax->GetValue().ToDouble(&new_value)) {
      bin_param1.rel_error_rho10_max = new_value; 
      }
      break;
    }
   case ID_BPARAM_ERTMAX:
    {
    if(TextCtrl_ertmax->GetValue().ToDouble(&new_value)) {
      bin_param1.error_theta10_deg_max = new_value; 
      }
      break;
    }
   case ID_BPARAM_RMIN:
    {
    if(TextCtrl_rmin->GetValue().ToDouble(&new_value)) {
      bin_param1.rho10_min = new_value; 
      }
      break;
    }
   case ID_BPARAM_NLRMIN:
    {
    if(TextCtrl_nlrmin->GetValue().ToDouble(&new_value)) {
      bin_param1.n_left_over_right_min = new_value; 
      }
      break;
    }
   case ID_BPARAM_NLRMAX:
    {
    if(TextCtrl_nlrmax->GetValue().ToDouble(&new_value)) {
      bin_param1.n_left_over_right_max = new_value; 
      }
      break;
    }
   case ID_BPARAM_NEGMAX:
    {
    if(TextCtrl_negmax->GetValue().ToDouble(&new_value)) {
      bin_param1.neg_max_for_patch = new_value; 
      }
      break;
    }
  }  // EOF switch
return;
}
