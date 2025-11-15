;
; Test julian()
; Arguments (inputs): 
;
month     = 1       ; <integer> month
day       = 9       ; <integer> day
year      = 1974    ; <integer> year
hour      = 9       ; <integer> hour
minute    = 0       ; <integer> minute
second    = 0       ; <integer> second
;
; Output variables must be mentioned before 
; the function call in order to reserve memory
; Arguments (outputs):
;
ierr      = 0       ; <integer> error flag (0 = OK)
date      = 0.0     ; <real> julian date             

print, julian(day,month,year,hour,minute,second,ierr,date)
;
if(ierr eq 0) then print, 'month   =', month 
if(ierr eq 0) then print, 'day     =', day
if(ierr eq 0) then print, 'year    =', year
if(ierr eq 0) then print, 'hour    =', hour
if(ierr eq 0) then print, 'minute  =', minute
if(ierr eq 0) then print, 'second  =', second
if(ierr eq 0) then print, 'This is julian date =', date
if(ierr ne 0) then print, 'julian ERROR :', ierr
;
; Test call_mcd()
; Arguments (inputs):
;
;TIME
datekey = 1       ; <integer> type of input date (1=Mars date)
ls = 120.0
localtime = 2.0   ; <real> local time (in martian hr) at lon. xlon
xdate = ls        ; <double precision> date (IF datekey = 1 : Value of Ls [deg.])
; LOCATION
zkey = 3          ; <integer> type of vertical coordinate xz (3 = above surface [m])
xz = 100.0        ; <real> vertical coordinate (m or Pa, depends on zkey)
; set hires flag
hireskey = 1      ; <integer> (1 = switch to high res. topography)
xlat = 5.0        ; <real> latitude (degrees north)
xlon = 6.0        ; <real> longitude (degrees east)
; DUST scenario
dset ='/path/to/MCD_DATA/'  ; <character*50> data set
scena = 1         ; <integer> scenario (1 = Climatology ave solar)
perturkey = 1     ; <integer>  perturbation type (1= none)
seedin   = 7.0    ; <real>
gwlength = 0.0    ; <real>  for small scale (ie: gravity wave) perturbations;
extvarkeys = LONARR(101) ; add 1 element because indexes in IDL start at 0
for i=1, 100 do extvarkeys[i] = 0 ; <integer> array output type (extvar(i) = 0 : don't compute)
;
; Output variables must be mentioned before 
; the function call in order to reserve memory
; Arguments (outputs):
;
meanvar = FLTARR(6) ; add 1 element because indexes in IDL start at 0
for i=1, 5 do meanvar[i] = 0.0 ; <real> mean unperturbed values (array of 5)
extvar = FLTARR(101) ; add 1 element because indexes in IDL start at 0
for i=1, 100 do extvar[i] = 0.0 ; <real>  extra variables (array of 100)
pres = 0.0        ; <real> atmospheric pressure (Pa)
dens = 0.0        ; <real> atmospheric density (kg/m^3)
temp = 0.0        ; <real> atmospheric temperature (K)
zonwind = 0.0     ; <real> zonal wind component (East-West)
merwind = 0.0     ; <real> meridional wind component (North-South)
seedout = 0.0     ; <real> current value of the seed of the random number generator
ierr = 0          ; <integer> error flag (0 = OK)
;
;
;
print, call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime, $
         dset,scena,perturkey,seedin,gwlength,extvarkeys, $
         pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ierr)
;
print, 'ierr   =', ierr
if(ierr eq 0) then print, 'pres   =', pres, 'Pa' 
if(ierr eq 0) then print, 'dens   =', dens, 'kg/m3'
if(ierr eq 0) then print, 'temp   =', temp, 'K'
if(ierr eq 0) then print, 'zonwind=', zonwind, 'm/s'
if(ierr eq 0) then print, 'merwind=', merwind, 'm/s'
if(ierr eq 0) then print, 'meanvar=', meanvar[1:5]
if(ierr eq 0) then print, 'extvar=', extvar[1:85]
if(ierr ne 0) then print, 'CALL_MCD ERROR :', ierr
