/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* jlp_spmeas_utils.cpp 
* To compute the theta zero for speckle interferometry
* using images of stars which moves because of the Earth rotation 
*
* Assumes that the stars move roughtly along the X axis
*
* JLP 
* Version 23-03-2020
-------------------------------------------------------------------*/
#include "jlp_spmeas.h"
#include "jlp_spmeas_utils.h"
#include "jlp_spmeas_latex.h"
#include "jlp_fitsio.h"         // JLP_RDFITS_2D_dble
#include "jlp_numeric.h"        // jlp_sort_array ...
#include "jlp_patch_set1.h"     // POLY_CIRC_PATCH, CREATE_NOISE...
#include "jlp_wx_ipanel_speckle.h" // speckle_patch_statistics ... 

#define DEBUG

/* defined here:
int JLP_LoadFITSImage(char *filename1, char *comments1,
                      double **dble_image1, int *nx1, int *ny1);
int JLP_SaveFITSImage(char *filename1, char *comments1,
                      double *dble_image1, int nx1, int ny1);
int JLP_FilterWithCircProfileToMaskedImage(double *dble_image1, 
                                           int nx1, int ny1, int n_for_patches);
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
                                          int ny1, char *latex_fname, 
                                          int *n_meas);
int JLP_ProcessSpeckleBinaryMeasurement(double *dble_image1, int nx1,
                                        int ny1, char *latex_fname, 
                                        int *n_meas);
int JLP_CosmeticPatch2(double *dble_image1, int nx1, int ny1,
                       double xc, double yc, double radius);
*/

/******************************************************************
* Load image from FITS file
*
* INPUT:
* filename1: name of the file to be read 
*
* OUTPUT:
* comments1: comments of the file that was read
* dble_image1: array with the data contained in FITS file
* nx1, ny1: size of the image 
*
*******************************************************************/
int JLP_LoadFITSImage(char *filename1, char *comments1,
                      double **dble_image1, int *nx1, int *ny1)
{
char errmess1[200];
int status, nz1, iframe;

// iframe: image plane to be loaded (from 1 to nz1)
iframe = 1;
status = JLP_RDFITS_2D_dble(dble_image1, nx1, ny1, &nz1, iframe, filename1,
                             comments1, errmess1);
if (status) {
  fprintf(stderr, "Couldn't load image from %s \n %s\n", filename1, errmess1);
  }

return(status);
}
/******************************************************************
* Save image to FITS files
*
* INPUT:
* filename1: name of the file to be read
*
* OUTPUT:
* comments1: comments of the file that was read
* dble_image1: array with the data contained in FITS file
* nx1, ny1: size of the image
*
*******************************************************************/
int JLP_SaveFITSImage(char *filename1, char *comments1,
                      double *dble_image1, int nx1, int ny1)
{
char errmess1[200];
int status;

status = JLP_WRFITS_2D_dble(dble_image1, nx1, ny1, nx1, filename1,
                             comments1, errmess1);
if (status) {
  fprintf(stderr, "Couldn't save image to %s \n %s\n", filename1, errmess1);
  }

return(status);
}
/*********************************************************************
* Filter with the minimum circular profile
* and generate a masked image
*
* INPUT
*  n_for_patches: number of points above the threshold
*
**********************************************************************/
int JLP_FilterWithCircProfileToMaskedImage(double *dble_image1, 
		                            int nx1, int ny1, int n_for_patches)
{
double *tmp_ima, thresh0;
int i, j, npts, iw, status, ixc, iyc, irho2;
#ifdef DEBUG
double  min0, max0;
char tmp_fname[128], tmp_comments[80];
#endif

npts = nx1 * ny1;
tmp_ima = new double[npts];
for(i=0; i < npts; i++) tmp_ima[i] = dble_image1[i];

#ifdef DEBUG
min0 = 1.e+8;
max0 = -1.e+8;
for(i=0; i < npts; i++) {
   if(dble_image1[i] > max0) max0 = dble_image1[i];
   if(dble_image1[i] < min0) min0 = dble_image1[i];
   }
// printf("DEBUG/min=%.2f max=%.2f\n", min0, max0);
sprintf(tmp_fname, "00_tmp_%d.fits", n_for_patches);
strcpy(tmp_comments, 
       "JLP_FilterWithCircProfileToMaskedImage/JLP_SubtractMinCircProfile");
JLP_SaveFITSImage(tmp_fname, tmp_comments, dble_image1, nx1, ny1);
#endif

ixc = nx1/2;
iyc = ny1/2;
status = JLP_SubtractMinCircProfile(tmp_ima, nx1, ny1, ixc, iyc);
for(i=0; i < npts; i++) dble_image1[i] = tmp_ima[i];

// Sort the intensity values of the image
 jlp_sort_array_intensities(dble_image1, nx1, ny1, tmp_ima, npts);

// n_for_patches: number of points above the threshold
// Threshold=total-n_for_patches values
iw = npts - n_for_patches;
if(iw < 0) iw = npts - 10;
thresh0 = tmp_ima[iw];
while((thresh0 <= 0.) && (iw < npts -1)) thresh0 = tmp_ima[iw++];
printf("JLP_FilterWithCircProfileToMaskedImage/Threshold=%.2f (npts-iw=%d)\n", thresh0, npts - iw);
for(i=0; i < nx1; i++)
  for(j=0; j < ny1; j++)
    {
    irho2 = SQUARE(i - ixc) + SQUARE(j - iyc);
// Discard all points too close to the center too:
    if((dble_image1[i + j * nx1] < thresh0) || (irho2 <= 9))
        dble_image1[i + j * nx1] = 0.;
    }

#ifdef DEBUG
min0 = 1.e+8;
max0 = -1.e+8;
for(i=0; i < npts; i++) {
   if(dble_image1[i] > max0) max0 = dble_image1[i];
   if(dble_image1[i] < min0) min0 = dble_image1[i];
   }
// printf("DEBUG/min=%.2f max=%.2f\n", min0, max0);
sprintf(tmp_fname, "0000_tmp_%d.fits", n_for_patches);
strcpy(tmp_comments, 
       "JLP_FilterWithCircProfileToMaskedImage/JLP_SubtractMinCircProfile");
JLP_SaveFITSImage(tmp_fname, tmp_comments, dble_image1, nx1, ny1);
#endif

delete[] tmp_ima;
return(0);
}
/**************************************************************************
* Compute the minimum circular profile of an image
*
* INPUT:
* image0: input 2D array
* nx0, ny0: dimension of the array
*
* OUTPUT:
* value0: value corresponding to the median
**************************************************************************/
int JLP_SubtractMinCircProfile(double *image0, int nx0, int ny0,
                                     int ixc0, int iyc0)
{
double *profile, rho;
int i, j, npts, nprof, iprof;

npts = nx0 * ny0;

// To avoid round off error, I multiply by 4 (and add 0.5):
nprof = 4 * (int)sqrt((double)(SQUARE(nx0)+SQUARE(ny0)));
profile = new double[nprof];

for(i = 0; i < nprof; i++) {
   profile[i] = 1.e+9;
   }

// Compute the minimum profile
  for(i = 0; i < nx0; i++)
    {
    for(j = 0; j < ny0; j++)
      {
      rho = SQUARE(i - ixc0) + SQUARE(j - iyc0);
      if(rho > 0.) rho = sqrt(rho);
// To avoid round off error, I multiply by 4 and add 0.5:
       iprof = (int)(4. * rho + 0.5);
       if(image0[i + j * nx0] < profile[iprof])
            profile[iprof] = image0[i + j * nx0];
      }
    }

// Subtract the minimum profile
  for(i = 0; i < nx0; i++)
    {
    for(j = 0; j < ny0; j++)
      {
      rho = SQUARE(i - ixc0) + SQUARE(j - iyc0);
      if(rho > 0.) rho = sqrt(rho);
// To avoid round off error, I multiply by 4 and add 0.5:
       iprof = (int)(4. * rho + 0.5);
       image0[i + j * nx0] -= profile[iprof];
      }
    }

delete[] profile;
return(0);
}
/*************************************************************************
* Measurement of binaries from circular profile (for autocorrelations)
**************************************************************************/
int JLP_BinaryMeasFromCircProfile(double *dble_image1, int nx1, int ny1,
                              int n_for_patches,
                              double *xc_left, double *yc_left, 
                              double *xc_right, double *yc_right,
                              double *radius,
                              double *rho10, double *theta10,
                              double *error_rho10, double *error_theta10,
                              double *n_left_over_right)
{
int status, ixc, iyc;

ixc = nx1/2;
iyc = ny1/2;
*xc_left = (double)ixc;
*yc_left = (double)iyc;
*xc_right = (double)ixc;
*yc_right = (double)iyc;
*radius = 10;
*rho10 = 0.;
*theta10 = 0.;
*error_rho10 = 0.;
*error_theta10 = 0.;
*n_left_over_right = 0.;

// Apply circular profile filter
JLP_FilterWithCircProfileToMaskedImage(dble_image1, nx1, ny1, 
		                           n_for_patches);

// Measure symmetric spots if any
status = JLP_MeasureBinaryFromMaskedImage(dble_image1, nx1, ny1, 
                                    rho10, theta10, error_rho10, error_theta10);

printf("Measurements: rho0=%.2f+/-%.2f pix. drho/rho=%.2f  theta0=%.2f+/-%.2f deg.\n",
       *rho10, *error_rho10, (*error_rho10)/(*rho10), *theta10 * (180. / PI),
       *error_theta10 * (180. / PI));

// err_rho/rho<0.15 too small, error_theta<15deg  too small
// Mini rho10 is 5 pixels to avoid too many false detections
// The scale is 0.0738 "/pixel, the diffraction limit is 0.16" or 2.16 pixels
if((*error_rho10 > 0.01) && (*error_theta10 * (180. / PI) > 0.01)
 && ((*error_rho10)/(*rho10) < 0.3) && (*error_theta10 * (180. / PI) < 30.)
 && (*rho10 > 5.)) {
      status = 0;
      } else {
      status = -1;
      }

// Determine the centers and the symmetry of the spots
if(status == 0) {
  status = JLP_XcentersFromMaskedImage(dble_image1, nx1, ny1, 
                                       *rho10, *theta10, xc_left, yc_left,
                                       xc_right, yc_right, n_left_over_right);
  printf("Xcenters: xc_right=%.1f yc_right=%.1f xc_left=%.1f yc_left=%.1f nleft_over_right=%.2f\n", 
      *xc_right, *yc_right, *xc_left, *yc_left, *n_left_over_right);
  if((*n_left_over_right > 0.2) && (*n_left_over_right < 5.)) {
      status = 0;
      } else {
      status = -1;
      }
  }

if(status == 0) {
// Radius increasing with distance to center:
// 10 at rho=70
  *radius = 3. + (*rho10) * 0.1;
  }

printf("Goodness of the fit: status=%d\n", status);

return(status);
}
/*************************************************************************
* Measurement of binaries from masked image 
*
* INPUT:
*
* OUTPUT:
*   rho10, err_rho10 (pixels)
*   theta10, error_theta10 (radians)
**************************************************************************/
int JLP_MeasureBinaryFromMaskedImage(double *dble_image1, int nx1, int ny1,
                               double *rho10, double *theta10,
                               double *error_rho10, double *error_theta10)
{
double min0, max0, rho2_min, rho2_max;
int i_iter, i, j, ixc, iyc;
double rho1, tan1, weight1, sum_weight1;
double r1, rr1, t1, tt1;

*rho10 = 0.;
*theta10 = 0.;
*error_rho10 = 0.;
*error_theta10 = 0.;

ixc = nx1 / 2;
iyc = ny1 / 2;

// Two iterations: first reject all points inside a circle of diameter=3
// then reject all points inside a circle of diameter rho/2
rho2_min = 9.;
rho2_max = SQUARE(nx1);
for( i_iter = 0; i_iter < 2; i_iter++) {

// Looking for the sky level as the minimum (far from the center):
min0 = 1.E+8;
max0 = -1.E+8;
for(i=0; i < nx1; i++)
  {
  for(j=0; j < ny1; j++)
    {
// Measure outside the circle centered on the center of the image, of radius=3
    if((SQUARE(i - ixc) + SQUARE(j - iyc) > rho2_min)
      && (SQUARE(i - ixc) + SQUARE(j - iyc) < rho2_max)
       && (dble_image1[i + j * nx1] > 0.)) {
       if (dble_image1[i + j * nx1] < min0) min0 = dble_image1[i + j * nx1];
       if (dble_image1[i + j * nx1] > max0) max0 = dble_image1[i + j * nx1];
       }
    }
  }
printf("JLP_MeasureBinaryFromMaskedImage/min=%.2f max=%.2f\n", min0, max0);

sum_weight1 = 0.;
r1 = 0.;
t1 = 0.;
rr1 = 0.;
tt1 = 0.;
for(i=0; i < nx1; i++)
  {
  for(j=0; j < ny1; j++)
    {
// Measure outside the circle centered on the center of the image, of radius=3
    if((SQUARE(i - ixc) + SQUARE(j - iyc) > rho2_min)
     && (SQUARE(i - ixc) + SQUARE(j - iyc) < rho2_max)
       && (dble_image1[i + j * nx1] > min0) && (i-ixc != 0.)) {
      rho1 = (double)(SQUARE(i - ixc) + SQUARE(j - iyc));
      tan1 = (double)(j - iyc)/(double)(i - ixc);
      if(rho1 > 0) {
        rho1 = sqrt(rho1);
        weight1 = dble_image1[i + j * nx1] - min0;
        r1 += rho1 * weight1;
        t1 += atan(tan1) * weight1;
        rr1 += SQUARE(rho1) * weight1;
        tt1 += SQUARE(atan(tan1)) * weight1;
        sum_weight1 += weight1;
        }
     }
  }
}

// Measurement and error:
if(sum_weight1 > 0.) {
 *rho10 = r1 / sum_weight1;
 *theta10 = t1 / sum_weight1;
 *error_rho10 = (rr1 / sum_weight1) - SQUARE(*rho10);
 *error_theta10 = (tt1 / sum_weight1) - SQUARE(*theta10);
 if(*error_rho10 > 0) *error_rho10 = sqrt(*error_rho10);
 if(*error_theta10 > 0) *error_theta10 = sqrt(*error_theta10);
}
/* DEBUG
*/
printf("JLP_MeasureBinaryFromMaskedImage/ rho0=%.2f+/-%.2f pix. theta0=%.2f+/-%.2f deg. (min0=%.3e sumweight=%.3e rho2_min=%.2f))\n",
   *rho10, *error_rho10, *theta10 * (180. / PI),
   *error_theta10 * (180. / PI), min0, sum_weight1, rho2_min);
// For the iteration i_iter=1, use a minimum radius of rho/3:
rho2_min = MINI(9., SQUARE(*rho10/3.));
rho2_max = SQUARE(*rho10 * 1.2);
}

return(0);
}
/*************************************************************************
* Measurement of binaries from masked image
* Determine the centers and the symmetry of the spots
*
* INPUT:
*   rho10  (pixels) 
*   theta10 (radians)
*
* OUTPUT:
**************************************************************************/
int JLP_XcentersFromMaskedImage(double *dble_image1, int nx1, int ny1,
                                double rho10, double theta10,
                                double *xc_left, double *yc_left,
                                double *xc_right, double *yc_right,
                                double *n_left_over_right)
{
double xcl, ycl, xcr, ycr, nl, nr; 
int i, j, ixc, iyc; 
double xvect, yvect, xvect0, yvect0, scalar_product;

ixc = nx1 / 2;
iyc = ny1 / 2;
xvect0 = rho10 * cos(theta10);
yvect0 = rho10 * sin(theta10);

nl = 0.;
nr = 0.;
xcl = 0.;
ycl = 0.;
xcr = 0.;
ycr = 0.;

for(i = 0; i < nx1; i++) {
  for(j = 0; j < ny1; j++) {
  if((dble_image1[i + j * nx1] > 0.)
     && (SQUARE(i - ixc) + SQUARE(j - iyc) > 9)) {
     xvect = (double)(i - ixc);
     yvect = (double)(j - iyc);
     scalar_product = xvect * xvect0 + yvect * yvect0;
     if(scalar_product > 0) {
       xcr += i - ixc;
       ycr += j - iyc;
       nr++;
     } else {
       xcl += i - ixc;
       ycl += j - iyc;
       nl++;
     }
   }
  } // EOF loop on j
} // EOF loop on i

*n_left_over_right = 0.;
if(nr > 0) {
  xcr /= nr;
  ycr /= nr; 
  *n_left_over_right = nl / nr;
}
if(nl > 0) {
  xcl /= nl;
  ycl /= nl; 
}
*xc_right = ixc + xcr;
*yc_right = iyc + ycr;
*xc_left = ixc + xcl;
*yc_left = iyc + ycl;

return(0);
}
/***********************************************************************
* Automatic speckle binary measurement
* Measure the position of the two autocorrelation peaks:
*  start with the left button is pressed by the user
*  and then measure the symmetric peak relative to the center of the frame
*
* INPUT:
*  latex_fname : latex file to be opened as append.
*  n_meas : number of measurements
***********************************************************************/
int JLP_AutomaticSpeckleBinaryMeasurement(double *original_image1, int nx1, 
                                          int ny1, char *fits_fname,
                                          char *latex_fname, 
                                          int *n_meas)
{
int bad_fit, status, output_poly_order;
int ifail, n_for_patches;
double *tmp_image1;
double mean_sky, sigma_sky, max_value;
double diam, xac, yac, flux, maxi;
char log_message[256], method[20], data_fname[128];
int i, i_iter, unresolved;
double xc_left, yc_left, xc_right, yc_right, radius, b_maxi, g_maxi;
double rho10, theta10, error_rho10, error_theta10, n_left_over_right;
FILE *fp_data;

strcpy(data_fname, "spmeas_data.tmp");
if((fp_data = fopen(data_fname,"w")) == NULL) {
  printf("JLP_AutomaticSpeckleBinaryMeasurement/Error writing data file >%s<\n",
          data_fname);
  return(-1);
  }

  tmp_image1 = new double[nx1 * ny1];
 
unresolved = 0;
for(i_iter = 0; i_iter < 2; i_iter++) {

// Apply process of each iteration to the original image:
 for(i = 0; i < nx1 * ny1; i++) tmp_image1[i] = original_image1[i];

 n_for_patches = 40 + i_iter * 10;
 printf("\niter=%d n_for_patches=%d\n", i_iter, n_for_patches);
// Circular profile filter:
 status = JLP_BinaryMeasFromCircProfile(tmp_image1, nx1, ny1, 
                                        n_for_patches, &xc_left, &yc_left, 
                                        &xc_right, &yc_right, &radius,
                                        &rho10, &theta10, &error_rho10, 
                                        &error_theta10, &n_left_over_right);
 if(status != 0) {
   unresolved = 1;
 } else {
   status = JLP_PatchProcessAndSaveBinaryMeasurements(fp_data, original_image1,
                                             nx1, ny1, xc_left, yc_left, 
                                             xc_right, yc_right, radius,
                                             rho10, theta10,
                                             error_rho10, error_theta10,
                                             n_meas);
   if(status != 0) {
    unresolved = 1; 
    break;
    }
 } // EOF unresolved==0

} // EOF loop on i_iter

// Close data file:
fclose(fp_data);

/* Write data file content to latex file: 
int JLP_speckle_process_data(char *data_fname,
                         char *latex_fname,
                         char *original_fits_fname,
                         char *processed_fits_fname, int unresolved);
*/
JLP_speckle_process_data(data_fname, latex_fname, fits_fname,
                         fits_fname, unresolved);

delete[] tmp_image1;
return(0);
}
/***********************************************************************
* Automatic speckle binary measurement
* Measure the position of the two autocorrelation peaks:
*  start with the left button is pressed by the user
*  and then measure the symmetric peak relative to the center of the frame
*
* INPUT:
*  latex_fname : latex file to be opened as append.
*  n_meas : number of measurements
***********************************************************************/
int JLP_PatchProcessAndSaveBinaryMeasurements(FILE *fp_data,
                                  double *original_image1, int nx1, int ny1,
                                  double xc_left, double yc_left, 
                                  double xc_right, double yc_right, 
                                  double radius,
                                  double rho10, double theta10,
                                  double error_rho10, double error_theta10,
                                  int *n_meas)
{
int bad_fit, status, output_poly_order;
int astrom_only, centered_polar_coordinates, ifail, n_for_patches;
double *tmp_image1;
double mean_sky, sigma_sky, max_value, negative_percent;
double diam, xc, yc, xac, yac, flux, maxi;
char log_message[256], method[20], data_fname[128];
int i, i_iter;
double b_maxi, g_maxi;

// astrom_only and centered_polar_coordinates set according to BinariesMode
  astrom_only = 1;
  centered_polar_coordinates = 1;

// Apply this to the original image:
tmp_image1 = new double[nx1 * ny1];

for(i_iter = 0; i_iter < 2; i_iter++) {

  if(i_iter == 0) {
    xc = xc_left;
    yc = yc_left;
  } else {
    xc = xc_right;
    yc = yc_right;
  }

for(i = 0; i < nx1 * ny1; i++) tmp_image1[i] = original_image1[i];

/**************************************************************
* xc, yc: center (in pixels) of the center of the circle
* radius: radius of the circle around the patch to be measured
* b_maxi: maximum in the circle (after subtraction of the sky background)
* g_maxi: maximum of the Gaussian (after subtraction of the sky background)
************************************************************************/
b_maxi = -1;
g_maxi = -1;
// Speckle: BinariesMode = 1
// Call JLP_Cosmetic Patch2 (with polynomial order = 3 by default
// and  m_polynomial_method = 0/1 (profile/polynomial)

 status = JLP_CosmeticPatch2(tmp_image1, nx1, ny1, xc, yc, radius);
/*** DEBUG
***/
printf("UUU OK: From JLP_CosmeticPatch2 with xc=%.1f yc=%.1f radius= %.1f status=%d\n", xc, yc, radius, status);

// Compute the difference (original image - new flattened image)
  for(i = 0; i < nx1 * ny1; i++) 
            tmp_image1[i] = original_image1[i] - tmp_image1[i];

/* Statistics on the edges: */
  speckle_patch_statistics(tmp_image1, nx1, ny1, xc, yc, radius, &mean_sky,
                           &sigma_sky, &max_value, &negative_percent, &bad_fit);
/*** DEBUG
*/
printf("UUU OK: From speckle_patch_statistics mean_sky=%.2f bad_fit=%d\n", 
        mean_sky, bad_fit);

 printf("Background: mean=%.2f, sigma=%.2f mean/max=%.2f\
 with %d%% of negative values\n",
mean_sky, sigma_sky, mean_sky/max_value, (int)negative_percent);
// Error if large oscillations in the background:
// 17 too small...
 if(negative_percent > 20.) {
    delete[] tmp_image1;
    return(-1);
  }

// Barycenter on the difference of the two images:
   diam = 2. * radius;
   status = astrom_barycenter(tmp_image1, nx1, ny1, xc, yc, diam, mean_sky,
                                &xac, &yac, &maxi, &flux);
/**** DEBUG
*/
printf("UUU OK: From astrom_barycenter xac=%.2f yac=%.2f status=%d\n", 
        xac, yac, status);

   if(status) {
     fprintf(stderr, "astrom_barycenter/Error: central patch is null!\
Empty circle or null sum\n");
   } else {
     b_maxi = maxi;
     strcpy(method,"Bary.");
     astrom_output_to_logfile(xac, yac, maxi, flux, xc, yc, diam,
                              output_poly_order, mean_sky, sigma_sky,
                              nx1, ny1, log_message, method,
                              astrom_only, centered_polar_coordinates);
printf("log_message: %s\n", log_message);
fprintf(fp_data, "%s\n", log_message);

// Increment the number of measurements:
   *n_meas++;
   } // EOF else

// /* Gaussian fit on the difference of the two images: */
// astrom_only = 1, centered_polar_coordinates=1
   diam = 2. * radius;
   status = astrom_gaussian_fit(tmp_image1, nx1, ny1, xc, yc, diam, mean_sky,
                                &xac, &yac, &maxi, &flux, &ifail);
   if(status) {
     sprintf(log_message,"astrom_gaussian_fit/ifail = %d \n", ifail);
   } else {
   g_maxi = maxi;
   strcpy(method,"Gauss.");
   astrom_output_to_logfile(xac, yac, maxi, flux, xc, yc, diam,
                            output_poly_order, mean_sky, sigma_sky,
                            nx1, ny1, log_message, method,
                            astrom_only, centered_polar_coordinates);
printf("log_message: %s\n", log_message);
fprintf(fp_data, "%s\n", log_message);

// Increment the number of measurements:
   *n_meas++;
   }  // EOF status == 0
} // EOF i_iter

delete[] tmp_image1;
return(0);
}
/***********************************************************************
* JLP_Cosmetic patch
*
* INPUT:
*  xc, yc: user coordinates of the center of the patch
*  radius: radius of the circular patch
*
* OUTPUT:
*   status = 1 if image was not modified
*          = 0 if image was modified
*          = -1 if problem in POLY_CIRC_PATCH
***********************************************************************/
int JLP_CosmeticPatch2(double *dble_image1, int nx1, int ny1,
                       double xc, double yc, double radius)
{
double *saved_data1;
double *work_image1, *noise_array, sigma_sky;
int noise_dim, Patch_Dlg_answer, status;
double m_radius_fact, m_sigma_noise;
int m_poly_order, m_polynomial_method;
char err_message[128];
register int i;

// Set default values for background estimation:
// m_poly_order = 3 here by default
m_radius_fact = 1.8;
m_sigma_noise = 0.3;
m_poly_order = 3;
m_polynomial_method = 1;

// Store initial image to "saved_data"
saved_data1 = new double[nx1 * ny1];
work_image1 = new double[nx1 * ny1];
for(i = 0; i < nx1 * ny1; i++) saved_data1[i] = dble_image1[i];

// Calling POLY_CIRC_PATCH to prepare a new version of the image:
/*
int POLY_CIRC_PATCH(double *image1, int nx1, int ny1, int idim1,
                    double xp, double yp, double diam, double diam_factor,
                    double *noise_array, int noise_dim, double sigma_noise,
                    int poly_order, double *sigma_sky);
int PROFILE_CIRC_PATCH(double *image1, int nx1, int ny1, int idim1,
                       double xp, double yp, double diam, double diam_factor,
                       double *noise_array, int noise_dim, double sigma_noise,
                       double *sigma_sky);
int CREATE_NOISE_ARRAY(double **noise_array, int noise_dim);
int DELETE_NOISE_ARRAY(double *noise_array);
*/
noise_dim = 256;

CREATE_NOISE_ARRAY(&noise_array, noise_dim);

// Retrieve original image:
for(i = 0; i < nx1 * ny1; i++) work_image1[i] = saved_data1[i];

// Process it:
if(m_polynomial_method) {
// Fit a polynomial to an annulus around the center
// and replace input values by noised computed values
  status = POLY_CIRC_PATCH(work_image1, nx1, ny1, nx1, xc, yc, 2.*radius,
                           m_radius_fact, noise_array, noise_dim,
                           m_sigma_noise, m_poly_order, &sigma_sky,
                           err_message);
  if(status) fprintf(stderr, "Error in Poly_circ_patch\n");
} else {
// Compute a mean profile (centered in nx1/2, ny1/2)
// and replace input values by noised computed values
// inside of the circle of diameter diam centered in (xp, yp)
  status = PROFILE_CIRC_PATCH(work_image1, nx1, ny1, nx1, xc, yc, 2.*radius,
                              m_radius_fact, noise_array, noise_dim,
                              m_sigma_noise, &sigma_sky, err_message);
  if(status) fprintf(stderr, "Error in Profile_circ_patch");
}

DELETE_NOISE_ARRAY(noise_array);

// copy back to dble_image1:
for(i = 0; i < nx1 * ny1; i++) dble_image1[i] = work_image1[i];

delete[] work_image1;
delete[] saved_data1;
// status = 1 if image was not modified
// status = 0 if image was modified
// status = -2 if problem in POLY_CIRC_PATCH
return(status);
}
