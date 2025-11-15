; call_mcd.pro
;
; This function performs checking of variables before a call
; is made to an external routine via CALL_EXTERNAL.
;
; The data types need to be the same as those expected by the
; external routine.
;
FUNCTION CALL_MCD, zkey1,xz1,xlon1,xlat1,hireskey1,datekey1,xdate1,localtime1, $
                   dset1,scena1,perturkey1,seedin1,gwlength1,extvarkeys1, $
                   pres1,dens1,temp1,zonwind1,merwind1,meanvar1,extvar1,seedout1,ier1

nextrvar   = N_ELEMENTS(extvarkeys1)
nmeanvar   = N_ELEMENTS(meanvar1)
extvarkeys = LONARR(nextrvar)
extvar     = FLTARR(nextrvar) 
meanvar    = FLTARR(nmeanvar)
;
; Inputs
;
zkey      = long(zkey1)
xz        = float(xz1)
xlon      = float(xlon1)
xlat      = float(xlat1)
hireskey  = long(hireskey1)
datekey   = long(datekey1)
xdate     = double(xdate1)
localtime = float(localtime1)
perturkey = long(perturkey1)
seedin    = float(seedin1)
gwlength  = float(gwlength1)
scena     = long(scena1)
extvarkeys= long(extvarkeys1)
;
; Insert str size in the first byte
str        = byte(dset1)
str_size   = N_ELEMENTS(str)
dset       = BYTARR(str_size+1)
FOR i=1,str_size DO dset[i]=str[i-1]
dset[0]    = byte(str_size)
;
; outputs
;
meanvar = float(meanvar1)
extvar  = float(extvar1)
seedout = float(seedout1)
pres    = float(pres1)
dens    = float(dens1)
temp    = float(temp1)
zonwind = float(zonwind1)
merwind = float(merwind1)
ier     = long(ier1)

result  = 1
;
; Call the external routine:
;
;HELP, zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,perturkey,seedin,gwlength,scena,extvarkeys,dset,meanvar,extvar,seedout,pres,dens,temp,zonwind,merwind,ier
result = CALL_EXTERNAL('call_mcd.so', 'call_mcd_wrapf_', $
         zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime, $
         dset,scena,perturkey,seedin,gwlength,extvarkeys, $
         pres,dens,temp,zonwind,merwind,meanvar,extvar,seedout,ier)

IF (result NE 1) THEN MESSAGE, 'Error calling external routine call_mcd'
;
; Return outputs:
;
meanvar1 = meanvar
extvar1  = extvar
seedout1 = seedout
pres1    = pres
dens1    = dens
temp1    = temp
zonwind1 = zonwind
merwind1 = merwind
ier1     = ier

RETURN, ier
END
