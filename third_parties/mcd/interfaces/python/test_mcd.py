#! /usr/bin/env python3

# This python program illustrates how the Mars Climate Database
# fortran routines can be called interactively from python

from fmcd import mcd
import numpy as np

# 1. Inputs:
dset = '/path/to/MCD_DATA/' # path to 'MCD_DATA'
perturkey = 1 # default to no perturbation
seedin = 0 # perturbation seed (unused if perturkey=1)
gwlength = 0 # Gravity Wave length for perturbations (unused if perturkey=1)

# 1.1 Dates
choice_date=input('Use Earth date (e) or Mars date (m)?')
if (choice_date == "e") :
  datekey=0 # Earth date
  loct=0 # local time must then also be set to zero
  day,month,year,hour,minute,second = \
  input("Enter date: day/month/year/hour/minute/second: ").split('/')
  day=int(day)
  month=int(month)
  year=int(year)
  hour=int(hour)
  minute=int(minute)
  second=int(second)
  # now call Julian routine to convert to julian date
  (ier,xdate)=mcd.julian(month,day,year,hour,minute,second)
  print(" Julian date %16.8f",xdate)
else :
  datekey=1 # Mars date
  xdate=float(input("Enter solar longitude Ls (deg.):"))
  loct=float(input("Local time (0 < time < 24)?"))

# 1.2 Vertical coordinate
zkey=int(input("Select verical coordinate type (1: distance to center of planet, 2: height above areoid, 3: height above surface, 4: Pressure) "))
if (zkey == 1) :
  xz=float(input("Enter distance to planet center (m) "))
if (zkey == 2) :
  xz=float(input("Enter altitude above areoid (m) "))
if (zkey == 3) :
  xz=float(input("Enter altitude above surface (m) "))
if (zkey == 4) :
  xz=float(input("Enter pressure value (Pa) "))

# high resolution mode
hrkey=int(input("High resolution? (1: yes, 0: no) "))

# 1.3 Position
lat = float(input('Latitude (deg)?'))
lon = float(input('Longitude (deg)?'))

# 1.4 Dust and solar scenario
print("Dust scenario?")
print("1= Climatology       typical Mars year dust scenario")
print("                     average solar EUV conditions")
print("2= Climatology       typical Mars year dust scenario")
print("                     minimum solar EUV conditions")
print("3= Climatology       typical Mars year dust scenario")
print("                     maximum solar EUV conditions")
print("4= dust storm        constant dust opacity = 5 (dark dust)")
print("                     minimum solar EUV conditions")
print("5= dust storm        constant dust opacity = 5 (dark dust)")
print("                     average solar EUV conditions")
print("6= dust storm        constant dust opacity = 5 (dark dust)")
print("                     maximum solar EUV conditions")
print("7= warm scenario     dustier than Climatology scenario")
print("                     maximum solar EUV conditions")
print("8= cold scenario     clearer than Climatology scenario")
print("                     minimum solar EUV conditions")
dust=int(input(''))

# 1.5 perturbations
perturkey=int(input("Perturbation? (1:none, 2: large scale, 3: small scale, 4: small+large, 5: n sigmas) "))
if (perturkey > 1) :
  seedin=int(input("seedin? (only matters if adding perturbations) "))
if ((perturkey == 3) or (perturkey == 4)) :
  gwlength=float(input("Gravity wave length? (for small scale perturbations) "))

# 1.6 extra outputs
# here we only implement an all-or-nothing case
extvarkey=int(input("Output the extra variables? (yes==1; no==0) "))
if (extvarkey == 0) :
  extvarkeys = np.zeros(100)
else :
  extvarkeys = np.ones(100)

# 2. Call MCD
(pres, dens, temp, zonwind, merwind, \
 meanvar, extvar, seedout, ierr) \
 = \
 mcd.call_mcd(zkey,xz,lon,lat,hrkey, \
 datekey,xdate,loct,dset,dust, \
 perturkey,seedin,gwlength,extvarkeys )

# 3. Write outputs
print("temperature is {:.0f} K, pressure is {:.0f} Pa, density is {:5.3e} kg/m3, zonal wind is {:.1f} m/s, meridional wind is {:.1f} m/s".format(temp,pres,dens,zonwind,merwind))

