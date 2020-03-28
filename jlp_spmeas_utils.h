/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* jlp_spmeas_utils.cpp 
* To compute the theta zero for speckle interferometry
* using images of stars which moves because of the Earth rotation 
*
* JLP 
* Version 23-03-2020
-------------------------------------------------------------------*/

int JLP_LoadFITSImage(char *filename1, char *comments1,
                      double **dble_image1, int *nx1, int *ny1);
int JLP_SaveFITSImage(char *filename1, char *comments1,
                      double *dble_image1, int nx1, int ny1);
int JLP_FilterWithCircProfileToMaskedImage(double *dble_image1, 
		                  int nx1, int ny1, int n_for_patches);
int JLP_SubtractMinCircProfile(double *image0, int nx0, int ny0,
                               int ixc0, int iyc0);
int JLP_BinaryMeasFromCircProfile(double *dble_image1, int nx1, int ny1,
                              int n_for_patches,
                              double *xc_left, double *yc_left, 
                              double *xc_right, double *yc_right,
                              double *radius,
                              double *rho10, double *theta10,
                              double *error_rho10, double *error_theta10,
                              double *n_left_over_right);
int JLP_PatchProcessAndSaveBinaryMeasurements(FILE *fp_data,
                                  double *original_image1, int nx1, int ny1,
                                  double xc_left, double yc_left, 
                                  double xc_right, double yc_right, 
                                  double radius,
                                  double rho10, double theta10,
                                  double error_rho10, double error_theta10,
                                  int *n_meas);
int JLP_MeasureBinaryFromMaskedImage(double *dble_image1, int nx1, int ny1,
                               double *rho10, double *theta10,
                               double *error_rho10, double *error_theta10);
int JLP_XcentersFromMaskedImage(double *dble_image1, int nx1, int ny1,
                                double rho10, double theta10,
                                double *xc_left, double *yc_left,
                                double *xc_right, double *yc_right,
                                double *n_left_over_right);
int JLP_AutomaticSpeckleBinaryMeasurement(double *dble_image1, int nx1,
                                          int ny1, char *fits_fname,
                                          char *latex_fname, int *n_meas);
int JLP_CosmeticPatch2(double *dble_image11, int nx11, int ny11,
                       double xc, double yc, double radius);
