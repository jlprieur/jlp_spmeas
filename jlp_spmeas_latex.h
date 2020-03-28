/********************************************************************
* jlp_spmeas_latex.h
* to save the data to LaTeX file
*
* JLP
* Version 28/03/2020
********************************************************************/
int JLP_speckle_process_data(char *data_fname,
                         char *latex_fname,
                         char *original_fits_fname,
                         char *processed_fits_fname,
                         int unresolved);
int read_int_from_cstring(char *cstring, char *mykey, int mykey_length,
                          int *ivalue, int *found);
int decode_info_from_FITS_file(const char *original_fits_fname, char *date, 
                               char *filter, char *object_name, double *epoch, 
                               double *year, int *xbin, int *ybin, 
                               char *comments);
