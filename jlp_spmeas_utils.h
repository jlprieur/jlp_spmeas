/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* jlp_spmeas_utils.cpp 
* To compute the theta zero for speckle interferometry
* using images of stars which moves because of the Earth rotation 
*
* JLP 
* Version 23-03-2020
-------------------------------------------------------------------*/
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>                   /* isprint... */

#include "jlp_spmeas_def.h"  // PI, MAXI, SQUARE, BIN_PARAM...

int process_list_for_jlp_spmeas(char *list_fname, char *latex_fname,
                                       char *lofile_csv);
int jlp_spmeas_for_one_file(char *in_name, char *latex_fname,
                                   char *logfile_csv);

int JLP_FilterWithCircProfileToMaskedImage(double *dble_image1, 
		                  int nx1, int ny1, int n_for_patches);
int JLP_SubtractMinCircProfile(double *image0, int nx0, int ny0,
                               int ixc0, int iyc0);
int JLP_BinaryMeasFromCircProfile(double *dble_image1, int nx1, int ny1,
                              int n_for_patches,
                              double *xc_left, double *yc_left, 
                              double *xc_right, double *yc_right,
                              double *rho10, double *theta10,
                              double *error_rho10, double *error_theta10,
                              double *n_left_over_right);
int JLP_BinaryMeasFirstCheck(double xc_left, double yc_left,
                             double xc_right, double yc_right,
                             double rho10, double theta10,
                             double error_rho10, double error_theta10,
                             double n_left_over_right, 
                             BIN_PARAM bin_param0, int *unresolved);
int JLP_Init_BIN_PARAM(BIN_PARAM *bin_param0);
int JLP_Load_BIN_PARAM_from_file(char *param_fname0, BIN_PARAM *bin_param0);
int JLP_Save_BIN_PARAM_to_file(char *param_fname0, BIN_PARAM bin_param0);
int JLP_PatchProcessAndSaveBinaryMeasurements(FILE *fp_data,
                                  double *original_image1, int nx1, int ny1,
                                  double xc_left, double yc_left, 
                                  double xc_right, double yc_right, 
                                  double radius0,
                                  double rho10, double theta10,
                                  double error_rho10, double error_theta10,
                                  double *negative_percent, int *n_meas,
                                  double *rho_mean, double *theta_mean);
int JLP_MeasureBinaryFromMaskedImage(double *dble_image1, int nx1, int ny1,
                               double *rho10, double *theta10,
                               double *error_rho10, double *error_theta10);
int JLP_XcentersFromMaskedImage(double *dble_image1, int nx1, int ny1,
                                double rho10, double theta10,
                                double *xc_left, double *yc_left,
                                double *xc_right, double *yc_right,
                                double *n_left_over_right);
int JLP_BinaryMeasInitCsvTable(char *logfile_csv);
int JLP_BinaryMeasToCsvTable(char *csv_fname, char *original_fits_name,
                             int n_for_patches, double rho_mean,
                             double theta_mean, double rho10, double theta10,
                             double error_rho10, double error_theta10,
                             double n_left_over_right, double negative_percent,
                             int unresolved, int *icheck, int ncheck);
int JLP_AutomaticSpeckleBinaryMeasurement(double *dble_image1, int nx1,
                                          int ny1, char *fits_fname,
                                          char *latex_fname, char *logfile_csv,
                                          BIN_PARAM bin_param0);
int JLP_CosmeticPatch2(double *dble_image11, int nx11, int ny11,
                       double xc, double yc, double radius0);
