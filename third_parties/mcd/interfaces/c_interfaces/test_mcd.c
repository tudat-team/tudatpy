#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

/**************************************************************************/
/* This program is an example to illustrate how the Mars Climate Database */
/* fortran routines can be called from a C program.                       */
/* Check the Makefile for compilation instructions and comments.          */
/**************************************************************************/
/* Checked to work with Gnu gcc compiler (version 4.1.2), as well as      */
/* with the Portland Group pgcc compiler (version 11.3) on Linux.         */
/* E. Millour, May 2006                                                   */
/* Adapted to MCD v4.2 "call_mcd" software, E. Millour 01/2007            */
/* Adapted to MCD v5.0 "call_mcd" software, E. Millour 10/2012            */
/* Adapted to MCD v6.1 "call_mcd" software, T. Pierron 06/2022            */
/**************************************************************************/
#include "mcd.h"

int main(){

char choice_date;  /* flag for date type */
int day,month,year,hour,min,sec; /* earth date */
float ls; /* solar longitude, for user input (if using mars date) */
int i;
/******* call_mcd arguments: *******/
/* call_mcd inputs: */
int zkey; /* flag to choose the type of z coordinates */
float xz; /* vertical coordinate value (m or Pa) */
float xlon; /* east longitude (degrees) */
float xlat; /* north latitude (degrees) */
int hireskey; /* high resolution flag (0: off, 1: on) */
int datekey; /* date flag (0: Earth date 1: Mars date) */
double xdate; /* Julian date (if datekey=0) or solar longitude(if datekey=1) */
float localtime; /* local time at longitude xlon (only if datekey=1) */
char dset[]="path/to/MCD_DATA"; /* path to datasets */
/* or dset[]=" "; to use default path "MCD_DATA/" (see call_mcd.F) */
int dust; /* dust and solar EUV scenario */
int perturkey; /* perturbation type */
float seedin; /* random generator seed and flag (if perturkey=1,2,3 or 4) */
              /* or coefficient to multiply std. dev. by (if perturkey=5) */
float gwlength; /* Gravity wave wavelength (needed if perturkey=3 or 4) */
int extvarkeys[100]; /* extra output variables (1: yes, 0: no) */
/* call_mcd outputs: */
float pres; /* atmospheric pressure */
float ro;   /* atmospheric density */
float temp; /* atmospheric temperature */
float u;    /* zonal wind */
float v;    /* meridional wind */
float meanvar[5]; /* unperturbed values of main meteorological variables */
float extvar[100]; /* extra output variables */
float seedout; /* current value of random generator seed index */
int ier; /* call_mcd_() returned status code (=0 if all went well) */


/* 1. Inputs */

/* 1.1. Time */
do {
  fprintf(stdout,"Use Earth date (e) or Mars date (m)?\n");
  fscanf(stdin,"%c",&choice_date);
} while ((choice_date!='e')&&(choice_date!='m'));

if (choice_date=='e') { /* Earth date */
  datekey=0;
  localtime=0; /*compulsary with earth dates */
  fprintf(stdout,"Enter date: day month year hour minute second\n");
  fscanf(stdin,"%d%d%d%d%d%d",&day,&month,&year,&hour,&min,&sec);
  fprintf(stdout," %2d/%2d/%4d %2d:%2d:%2d\n",day,month,year,hour,min,sec);
  /* convert to julian date */
  __mcd_MOD_julian(&month, &day, &year, &hour, &min, &sec, &ier, &xdate);
  fprintf(stdout,"Julian date: %16.8f\n",xdate);
}
else { /* Mars date */
  datekey=1;
  fprintf(stdout,"Enter solar longitude Ls (deg.):\n");
  fscanf(stdin,"%g",&ls);
  fprintf(stdout,"Local time (0 < time < 24)?\n");
  fscanf(stdin,"%g",&localtime);
  xdate=ls;
}

/* 1.2. Vertical coordinate */

do {
  fprintf(stdout,"Select the vertical coordinate type:\n");
  fprintf(stdout,"1: Distance to center of planet (meters)\n");
  fprintf(stdout,"2: height, above areoid (meters)\n");
  fprintf(stdout,"3: height, above surface (meters)\n");
  fprintf(stdout,"4: Pressure level (Pa)\n");
  fscanf(stdin,"%d",&zkey);

  if (zkey==1) {
    fprintf(stdout,"Enter distance to planet center (m)\n");
    fscanf(stdin,"%g",&xz);
  }
  if (zkey==2) {
    fprintf(stdout,"Enter above areoid height (m)\n");
    fscanf(stdin,"%g",&xz);
  }
  if (zkey==3) {
    fprintf(stdout,"Enter above surface height (m)\n");
    fscanf(stdin,"%g",&xz);
  }
  if (zkey==4) {
    fprintf(stdout,"Enter pressure value (Pa)\n");
    fscanf(stdin,"%g",&xz);
  }
} while ((zkey<1)||(zkey>4));

/* High resolution mode? */
do {
  fprintf(stdout,"High resolution? (1: yes, 0: no)\n");
  fscanf(stdin,"%d",&hireskey);
} while ((hireskey<0)||(hireskey>1));

/* 1.3. Position */

fprintf(stdout,"Latitude (deg)?\n");
fscanf(stdin,"%g",&xlat);
fprintf(stdout,"East longitude (deg)?\n");
fscanf(stdin,"%g",&xlon);

/* 1.4. Dust and solar scenario */

do {
fprintf(stdout,"Dust scenario?\n");
fprintf(stdout,"1= Climatology       typical Mars year dust scenario\n");
fprintf(stdout,"                     average solar EUV conditions\n");
fprintf(stdout,"2= Climatology       typical Mars year dust scenario\n");
fprintf(stdout,"                     minimum solar EUV conditions\n");
fprintf(stdout,"3= Climatology       typical Mars year dust scenario\n");
fprintf(stdout,"                     maximum solar EUV conditions\n");
fprintf(stdout,"4= dust storm        constant dust opacity = 4\n");
fprintf(stdout,"                     min solar EUV conditions\n");
fprintf(stdout,"5= dust storm        constant dust opacity = 4\n");
fprintf(stdout,"                     ave solar EUV conditions\n");
fprintf(stdout,"6= dust storm        constant dust opacity = 4\n");
fprintf(stdout,"                     max solar EUV conditions\n");
fprintf(stdout,"7= warm scenario     warm scenario: dustier than MY24\n");
fprintf(stdout,"                     max solar EUV conditions\n");
fprintf(stdout,"8= cold scenario     cold scenario: clearer than MY24\n");
fprintf(stdout,"                     min solar EUV conditions\n");
fscanf(stdin,"%d",&dust);
} while ((dust<1)||(dust>8));

/* 1.5 Perturbations */

fprintf(stdout,"Perturbation: none = 1 ;  large scale = 2 ;");
fprintf(stdout," small scale= 3 ; small+large = 4 ; n sig =5\n");
fscanf(stdin,"%d",&perturkey);
fprintf(stdout,"seedin and gwlength\n");
fscanf(stdin,"%g %g",&seedin,&gwlength);

/* 1.6 extra output */
do {
  /* here we only implement the case of either none or all output variables */
  fprintf(stdout,"Output the extra variables? (yes==1; no==0)\n");
  fscanf(stdin,"%d",&extvarkeys[0]);
  for (i=1;i<100;i++) extvarkeys[i]=extvarkeys[0] ;
}while ((extvarkeys[0]!=0)&&(extvarkeys[0]!=1));


/* 2. call fortran routine call_mcd */

__mcd_MOD_call_mcd(&zkey,&xz,&xlon,&xlat,&hireskey,
          &datekey,&xdate,&localtime,dset,&dust,
          &perturkey,&seedin,&gwlength,extvarkeys,
          &pres,&ro,&temp,&u,&v,meanvar,extvar,&seedout,&ier,
          strlen(dset));

/* 3. write output */

if(ier==0) {
  fprintf(stdout,"p   =             %g Pa\n",pres);
  fprintf(stdout,"rho =             %g kg/m^3\n",ro);
  fprintf(stdout,"T   =             %g K\n",temp);
  fprintf(stdout,"Zonal wind      = %g m/s\n",u);
  fprintf(stdout,"Meridional wind = %g m/s\n",v);

/* Note: In order to stick to the fortran-based description
   of call_mcd outputs (see the user's manual),
   printed indexes are shifted. */
  for (i=0;i<5;i++) {
    fprintf(stdout,"Meanvar (%d)=%g\n",i+1,meanvar[i]);
  }
  if (extvarkeys[0]!=0) { /* write all 76 extra variables */
    for (i=0;i<85;i++) {
      fprintf(stdout,"Extvar (%d)=%g\n",i+1,extvar[i]);
    }
  }
  else { /* write only the first 7 variables */
    for (i=0;i<13;i++) {
      fprintf(stdout,"Extvar (%d)=%g\n",i+1,extvar[i]);
    }
  }
}
else {
  fprintf(stdout,"CALL_MCD ERROR !!\n");
  fprintf(stdout,"         returned error code: %d\n",ier);
}

return 0;
}
