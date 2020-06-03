/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* jlp_spmeas_autom : routine for automatic measurements
*
* JLP 
* Version 13-03-2020
-------------------------------------------------------------------*/

#include "spm_frame.h"
#include "jlp_spmeas_utils.h"
#include "jlp_fitsio.h"

#ifdef MAIN
int main(int argc, char *argv[])
{
int status;
char list_fname[128], latex_fname[128], csv_fname[128];

printf(" Program jlp_spmeas to compute theta zero \n");
printf(" JLP Version 06-03-2020 \n");

/* One parameters only is allowed to run the program: */
/* Carefull: 7 parameters always, using JLP "runs" */
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 3)
  {
  printf(" Syntax: jlp_spmeas input_list output_latex_file\n"); 
  printf(" input_list: list of image fits files\n");
  printf(" Fatal: Syntax error: argc=%d\n",argc);
  exit(-1);
  }

/* Input of parameters with the command line: */
 strcpy(list_fname, argv[1]);
 strcpy(latex_fname, argv[2]);
 strcpy(csv_fname, "jlp_spmeas_csv.txt");

 status = process_list_for_jlp_spmeas(list_fname, latex_fname, csv_fname);

return(status);
}
#endif

/*************************************************************************
*
*************************************************************************/
int SpmFrame::ProcessListOfAutocFilesForAutom()
{
wxString str1;
int status, k, nfiles;
char in_fname[128];
FILE *fp_in;

if(list_fname1[0] == '\0') {
  wxMessageBox(_T("Error: please select a name for the input autoc list file"),
               _T("ProcessListOfAutocFilesForAutom"), wxICON_ERROR);
  return(-1);
  }


if(csv_fname1[0] == '\0') {
  wxMessageBox(_T("Warning: output csv file not set! I will set to jlp_spmeas_csv.txt"),
               _T("ProcessListOfAutocFilesForAutom"), wxICON_EXCLAMATION);
  strcpy(csv_fname1, "jlp_spmeas_csv.txt");
  }

if(latex_fname1[0] == '\0') {
  wxMessageBox(_T("Warning: output latex file not set! I will set to jlp_spmeas.tex"),
               _T("ProcessListOfAutocFilesForAutom"), wxICON_EXCLAMATION);
  strcpy(latex_fname1, "jlp_spmeas.tex");
  }

// Start by opening csv file:
JLP_BinaryMeasInitCsvTable(csv_fname1);

if((fp_in = fopen(list_fname1,"r")) == NULL) {
  str1.Printf(_T("Error reading list of files >%s<"), list_fname1);
  wxMessageBox(str1, _T("ProcessListOfAutocFilesForAutom"),
               wxICON_ERROR);
  return(-1);
  }

/* Loop on all the files: */
k = 0;
while(!feof(fp_in)) {
  if(fscanf(fp_in, "%s\n", in_fname)) {
   if((in_fname[0] != '%') && (in_fname[0] != '#')) {
     status  = AutoMeasureOfAutoc(in_fname);
     if(status == 0) k++;
     }
   }
}
fclose(fp_in);
nfiles = k;
// Save binary parameters:
  strcpy(param_fname1, "jlp_spmeas_param.txt");
 JLP_Save_BIN_PARAM_to_file(param_fname1, bin_param1);

str1.Printf(_T(" %d file(s) processed from input list (%s) \n Csv output: %s\n, LateX output: %s\nParameters saved to %s"), 
		nfiles, list_fname1, csv_fname1, latex_fname1, param_fname1);
wxMessageBox(str1, _T("ProcessListOfAutocFilesForAutom"), wxOK);

return(status);
}
/*************************************************************************
*
*************************************************************************/
int SpmFrame::AutoMeasureOfAutoc(char *in_name)
{
char comments1[128];
double *dble_image1;
int status, nx1, ny1;
wxString str1;

dble_image1 = NULL;

printf("\n*******************************************************\n");
printf("jlp_spmeas_for_one_file: processing %s\n", in_name);
status = JLP_LoadFITSImage(in_name, comments1, &dble_image1, &nx1, &ny1);

if(status != 0) {
    str1.Printf(_T("Error in JLP_LoadFITSImage loading %s, status=%d\n"), 
                in_name, status); 
    wxMessageBox(str1, _T("AutoMeasureOfAutoc"), wxICON_ERROR);
    return(-1);
    }

status = JLP_AutomaticSpeckleBinaryMeasurement(dble_image1, nx1, ny1, in_name, 
                                               latex_fname1, csv_fname1,
                                               bin_param1);

if(dble_image1 != NULL) delete[] dble_image1;
return(0);
}
