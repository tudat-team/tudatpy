#include <iostream>
#include <cmath>
#include <string>

////////////////////////////////////////////////////////////////////////////
// This program is an example to illustrate how the Mars Climate Database //
// fortran routines can be called from a C++ program.                     //
// Check the Makefile for compilation instructions and comments.          //
////////////////////////////////////////////////////////////////////////////
// Checked to work with Gnu g++ compiler (version 4.1.2), as well as      //
// with the Portland Group pgCC compiler (version 11.3) on Linux.          //
// Initial version freely and kindly distributed by John Underwood,       //
// Vorticity Ltd.                                                         //
// Ehouarn Millour, May 2006                                              //
// Adapted to MCD v4.2 "call_mcd" software, E. Millour 01/2007            //
// Adapted to MCD v5.0 "call_mcd" software, E. Millour 10/2012            //
////////////////////////////////////////////////////////////////////////////

#include "mcd.h"

using namespace std;

int main() {
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
string dset="/path/to/MCD_DATA/"; /* path to datasets */
/* or dset=" "; to use default path "MCD_DATA/" (see call_mcd.F) */
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


// 1. Inputs

// 1.1. Time
do {
  cout<<"Use Earth date (e) or Mars date (m)?"<<endl;
  cin>>choice_date;
} while ((choice_date!='e')&&(choice_date!='m'));

if (choice_date=='e') { /* Earth date */
  datekey=0;
  localtime=0; // compulsary with earth dates
  cout<<"Enter date: day month year hour minute second"<<endl;
  cin>>day>>month>>year>>hour>>min>>sec;
  cout<<day<<"/"<<month<<"/"<<year<<" "<<hour<<":"<<min<<":"<<sec<<endl;
  /* convert to julian date */
  __mcd_MOD_julian(&month, &day, &year, &hour, &min, &sec, &ier, &xdate);
  cout.setf(ios::fixed);
  cout<<"Julian date: "<<xdate<<endl;
  cout.unsetf(ios::fixed);
}
else { /* Mars date */
  datekey=1;
  cout<<"Enter solar longitude Ls (deg.):"<<endl;
  cin>>ls;
  cout<<"Local time (0 < time < 24)?"<<endl;
  cin>>localtime;
  xdate=ls;
}

// 1.2. Vertical coordinate
do {
  cout<<"Select the vertical coordinate type:"<<endl;
  cout<<"1: Distance to center of planet (meters)"<<endl;
  cout<<"2: height, above areoid (meters)"<<endl;
  cout<<"3: height, above surface (meters)"<<endl;
  cout<<"4: Pressure level (Pa)"<<endl;
  cin>>zkey;

  if (zkey==1) {
    cout<<"Enter distance to planet center (m)"<<endl;
    cin>>xz;
  }
  if (zkey==2) {
    cout<<"Enter above areoid height (m)"<<endl;
    cin>>xz;
  }
  if (zkey==3) {
    cout<<"Enter above surface height (m)"<<endl;
    cin>>xz;
  }
  if (zkey==4) {
    cout<<"Enter pressure value (Pa)"<<endl;
    cin>>xz;
  }
} while ((zkey<1)||(zkey>4));

// High resolution mode?
do {
  cout<<"High resolution? (1: yes, 0: no)"<<endl;
  cin>>hireskey;
} while ((hireskey<0)||(hireskey>1));

// 1.3. Position
cout<<"Latitude (deg)?"<<endl;
cin>>xlat;
cout<<"East longitude (deg)?"<<endl;
cin>>xlon;

// 1.4. Dust and solar scenario

do {
cout<<"Dust scenario?"<<endl;
cout<<"1= Climatology       typical Mars year dust scenario"<<endl;
cout<<"                     average solar EUV conditions"<<endl;
cout<<"2= Climatology       typical Mars year dust scenario"<<endl;
cout<<"                     minimum solar EUV conditions"<<endl;
cout<<"3= Climatology       typical Mars year dust scenario"<<endl;
cout<<"                     maximum solar EUV conditions"<<endl;
cout<<"4= dust storm        constant dust opacity = 4"<<endl;
cout<<"                     min solar EUV conditions"<<endl;
cout<<"5= dust storm        constant dust opacity = 4"<<endl;
cout<<"                     ave solar EUV conditions"<<endl;
cout<<"6= dust storm        constant dust opacity = 4"<<endl;
cout<<"                     max solar EUV conditions"<<endl;
cout<<"7= warm scenario     warm scenario: dustier than MY24"<<endl;
cout<<"                     max solar EUV conditions"<<endl;
cout<<"8= cold scenario     cold scenario: clearer than MY24"<<endl;
cout<<"                     min solar EUV conditions"<<endl;
cin>>dust;
} while ((dust<1)||(dust>8));

// 1.5 Perturbations
cout<<"Perturbation: none = 1 ;  large scale = 2 ;";
cout<<" small scale= 3 ; small+large = 4 ; n sig =5"<<endl;
cin>>perturkey;
cout<<"seedin and gwlength"<<endl;
cin>>seedin>>gwlength;

// 1.6 extra output 
do {
  // here we only implement the case of either none or all output variables
  cout<<"Output the extra variables? (yes==1; no==0)"<<endl;
  cin>>extvarkeys[0];
  for (i=0;i<100;i++) extvarkeys[i]=extvarkeys[0] ;
}while ((extvarkeys[0]!=0)&&(extvarkeys[0]!=1));


// 2. call fortran routine call_mcd 
__mcd_MOD_call_mcd(&zkey,&xz,&xlon,&xlat,&hireskey,
          &datekey,&xdate,&localtime,dset.c_str(),&dust,
          &perturkey,&seedin,&gwlength,extvarkeys,
          &pres,&ro,&temp,&u,&v,meanvar,extvar,&seedout,&ier,
          dset.length());

// 3. write output

if(ier==0) {
  cout<<"p   =             "<<pres<<" Pa"<<endl;
  cout<<"rho =             "<<ro<<" kg/m^3"<<endl;
  cout<<"T   =             "<<temp<<" K"<<endl;
  cout<<"Zonal wind      = "<<u<<" m/s"<<endl;
  cout<<"Meridional wind = "<<v<<" m/s"<<endl;

// Note: In order to stick to the fortran-based description
// of atmemcd outputs (see the user's manual),
// printed indexes are shifted. 
  for (i=0;i<5;i++) {
    cout<<"Meanvar ("<<i+1<<")="<<meanvar[i]<<endl;
  }
  if(extvarkeys[0]!=0) { /* write all 76 extra variables */
    for (i=0;i<85;i++) {
      cout<<"Extvar ("<<i+1<<")="<<extvar[i]<<endl;
    }
  }
  else { /* write only the first 7 variables */
    for (i=0;i<13;i++) {
      cout<<"Extvar ("<<i+1<<")="<<extvar[i]<<endl;
    }
  }
}
else {
  cout<<"CALL_MCD ERROR !!"<<endl;
  cout<<"         returned error code: "<<ier<<endl;
}

return 0;
}
