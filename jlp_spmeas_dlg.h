/******************************************************************************
* jlp_spmeas_dlg.h
* Dialog box used for selecting the binary parameters for automatic meaurements
*
* Author:      JLP 
* Version:     18/05/2020
******************************************************************************/
#ifndef jlp_spmeas_dlg_h    // sentry 
#define jlp_spmeas_dlg_h

#include "jlp_spmeas_def.h"  // BIN_PARAM

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wx/image.h"
#include "wx/file.h"
#include "wx/filename.h"
#include "wx/mstream.h"
#include "wx/wfstream.h"
#include "wx/quantize.h"
#include "wx/dcmemory.h"         // wxMemoryDC

/********************************************************************
* Class JLP_Spmeas_Dlg
*********************************************************************/

class JLP_Spmeas_Dlg: public wxDialog
{
public:

// Constructor:
     JLP_Spmeas_Dlg(wxFrame *parent, BIN_PARAM bin_param0, 
                    wxString title);
// Destructor: 
    ~JLP_Spmeas_Dlg(){
       };

    void InitPanel();

// Handling events:
     void OnOKButton( wxCommandEvent &event );
     void OnCancelButton( wxCommandEvent &event );
     void OnDefaultButton( wxCommandEvent &event );
     void OnChangeParam( wxCommandEvent& event );

// Accessors:
   int RetrieveData(BIN_PARAM *bin_param0) {
       bin_param0->error_theta10_deg_max = bin_param1.error_theta10_deg_max; 
       bin_param0->rel_error_rho10_max = bin_param1.rel_error_rho10_max; 
       bin_param0->rho10_min = bin_param1.rho10_min; 
       bin_param0->n_left_over_right_min = bin_param1.n_left_over_right_min; 
       bin_param0->n_left_over_right_max = bin_param1.n_left_over_right_max; 
       bin_param0->neg_max_for_patch = bin_param1.neg_max_for_patch; 
       return(0);
       }

protected:
   int DataIsOK() {
       if(bin_param1.rel_error_rho10_max < 0. 
          || bin_param1.error_theta10_deg_max < 0.) return(false); 
       if(bin_param1.rho10_min < 0. 
          || bin_param1.n_left_over_right_min < 0.
          || bin_param1.n_left_over_right_max < 0.
          || bin_param1.neg_max_for_patch < 0) return(false); 
       if(bin_param1.n_left_over_right_min > bin_param1.n_left_over_right_max) 
         return(false); 
       return(true);
       }

private:
    int initialized;
    BIN_PARAM bin_param1;
    wxTextCtrl *TextCtrl_ertmax; 
    wxTextCtrl *TextCtrl_relerrmax, *TextCtrl_rmin, *TextCtrl_nlrmin;
    wxTextCtrl *TextCtrl_nlrmax, *TextCtrl_negmax;

    DECLARE_EVENT_TABLE()
};

#endif               // EOF sentry
