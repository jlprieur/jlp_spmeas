/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Program jlp_spmeas 
* To compute the theta zero for speckle interferometry
* using images of stars which moves because of the Earth rotation 
*
* Assumes that the stars move roughtly along the X axis
*
* JLP 
* Version 13-03-2020
-------------------------------------------------------------------*/
#define DEBUG

#include "jlp_spmeas.h"
#include "jlp_spmeas_utils.h"


#define NCOEFF 2

/* Defined here: */
static int process_list_for_jlp_spmeas(char *list_fname, char *latex_fname);
static int jlp_spmeas_for_one_file(char *in_name, char *latex_fname);

int main(int argc, char *argv[])
{
int status;
char list_fname[128], latex_fname[128];

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

 status = process_list_for_jlp_spmeas(list_fname, latex_fname);

return(status);
}
/*************************************************************************
*
*************************************************************************/
static int process_list_for_jlp_spmeas(char *list_fname, char *latex_fname)
{
int status, k, nfiles;
char in_fname[128];
FILE *fp_in, *fp_out;

if((fp_in = fopen(list_fname,"r")) == NULL) {
  printf("process_list_for_jlp_spmeas/Error reading list of files >%s<\n",
          list_fname);
  return(-1);
  }

/* First go to count the files: */
k = 0;
while(!feof(fp_in)) {
  if(fscanf(fp_in, "%s\n", in_fname)) {
   if((in_fname[0] != '%') && (in_fname[0] != '#')) {
     k++;
     status  = jlp_spmeas_for_one_file(in_fname, latex_fname);
     }
   }
}
fclose(fp_in);
nfiles = k;
printf("%d files have been processed in this input list (%s)\n", k, list_fname);

return(status);
}
/*************************************************************************
*
*************************************************************************/
static int jlp_spmeas_for_one_file(char *in_name, char *latex_fname) 
{
char comments1[128];
double *dble_image1;
int status, nx1, ny1, n_meas;

printf("\n*******************************************************\n");
printf("jlp_spmeas_for_one_file: processing %s\n", in_name);
status = JLP_LoadFITSImage(in_name, comments1, &dble_image1, &nx1, &ny1);

if(status != 0) {
    printf("jlp_spmeas_for_one_file/Error in JLP_LoadFITSImage loading %s, status=%d\n", 
          in_name, status); 
    return(-1);
    }

n_meas = 0;
status = JLP_AutomaticSpeckleBinaryMeasurement(dble_image1, nx1, ny1, in_name, 
                                               latex_fname, &n_meas);

if(dble_image1 != NULL) delete[] dble_image1;
return(0);
}
