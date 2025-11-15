function [ires] = test_mlb2mcd()
%
% mlb_julian()
% Arguments (inputs): 
%
month     = 1;       % <integer> month
day       = 9;       % <integer> day
year      = 1974;    % <integer> year
hour      = 9;       % <integer> hour
minute    = 0;       % <integer> minute
second    = 0;       % <integer> second
%
% Output variables must be mentioned before 
% the function call in order to reserve memory
% Arguments (outputs):
%
for i=1:2
 ierr(i)      = 0;       % <integer> error flag (0 = OK)
 date(i)      = 0.0;     % <real> julian date
end             

mlb_julian(month,day,year,hour,minute,second,ierr,date);

if ierr(1)==0
 fprintf('month  = %g\n', month );
 fprintf('day    = %g\n', day   );
 fprintf('year   = %g\n', year  );
 fprintf('hour   = %g\n', hour  );
 fprintf('minute = %g\n', minute);
 fprintf('second = %g\n', second);
 fprintf('This is julian date = %g\n\n', date(1));
else
 fprintf('julian ERROR !!\n');
 fprintf('returned error code: %d\n\n',ierr(1));
end

%
% mlb_call_mcd()
% Arguments (inputs):
%
%TIME
datekey   = 1;       % <integer> type of input date (1=Mars date)
ls        = 120.0;   
localtime = 2.0;     % <real> local time (in martian hr) at lon. xlon
xdate     = ls;      % <double precision> date (IF datekey = 1 : Value of Ls [deg.])
% LOCATION
zkey      = 3;       % <integer> type of vertical coordinate xz (3 = above surface [m])
xz        = 100.0;   % <real> vertical coordinate (m or Pa, depends on zkey)
% set hires flag
hireskey  = 1;       % <integer> (1 = switch to high res. topography)
xlat      = 5.0;     % <real> latitude (degrees north)
xlon      = 6.0;     % <real> longitude (degrees east)

% DUST scenario
dset      = '/path/to/mcd/data/'; % <character*50> data set
dust      = 1;       % <integer> scenario (1 = Climatology ave solar)
perturkey = 1;       % <integer>  perturbation type (1= none)
seedin    = 7.0;     % <real> 
gwlength  = 0.0;     % <real>  for small scale (ie: gravity wave) perturbations;
extvarkeys(1)=0;     % <integer> array output type (extvar(i) = 0 : don't compute)

for i=2:100
 extvarkeys(i)=extvarkeys(1);
end % i
%
% Output variables must be mentioned before 
% the function call in order to reserve memory
% Arguments (outputs):
%
for i=1:5
 meanvar(i) = 0.0; % <real> mean unperturbed values (array of 5)
end

for i=1:100
 extvar(i) = 0.0; % <real>  extra variables (array of 100)
end

for i=1:2
 pres(i)    = 0.0; % <real> atmospheric pressure (Pa)
 ro(i)      = 0.0; % <real> atmospheric density (kg/m^3)
 temp(i)    = 0.0; % <real> atmospheric temperature (K)
 u(i)       = 0.0; % <real> zonal wind component (East-West)
 v(i)       = 0.0; % <real> meridional wind component (North-South)
 seedout(i) = 0.0; % <real> current value of the seed of the random number generator
 ier(i)     = 0;   % <integer> error flag (0 = OK)
end
%
%
%
mlb_call_mcd(zkey,xz,xlon,xlat,hireskey,datekey,xdate,localtime,dset,dust,perturkey,seedin,gwlength,extvarkeys,pres,ro,temp,u,v,meanvar,extvar,seedout,ier);

if ier(1)==0
 fprintf('p   = %g Pa\n', pres(1));
 fprintf('rho = %g kg/m3\n',ro(1));
 fprintf('T   = %g K \n',temp(1));
 fprintf('Zonal wind      = %g m/s\n',u(1));
 fprintf('Meridional wind = %g m/s\n',v(1));
 fprintf('\n');
 for i=1:5
  fprintf('meanvar(%d)= %g\n',i,meanvar(i));
 end
 fprintf('\n');         
 if extvarkeys(1)~=0 
  for i = 1:85
   fprintf('extvar(%d) = %g\n',i,extvar(i));
  end
 else
  for i=1:7
   fprintf('extvar(%d) = %g\n',i,extvar(i));
  end
 end
else
 fprintf('CALL_MCD ERROR !!\n');
 fprintf('returned error code: %d\n',ier(1));
end

ires = 0;

