clear all
format long
%
degrad = pi/180;
raddeg = 180/pi;
%
Re     = 6.378137e6; % Equatorial radius Earth (m)
omcb   = 7.292115e-5; % Rotational rate Earth (rad/s)
THETA  = 0.0; % I-frame is R-frame at t = 0 (Greenwhich Hour Angle)
%%
% Input relative spherical state (position: r, tau (lon), delta (lat))
%
h      =  119961.097;
r      =  Re + h;
tau    =  0.0;%-105.965657*degrad;
delta  =  90.0*degrad;%-22.0604427*degrad;
% Velocity (modulus Vrel, flight-path angle gamma, heading chi)
Vrel   =  7438.14752;
gamma  =  0.0*-1.42632303*degrad;
chi    =  0.0*70.4369784*degrad;
%
% Transformtion V-frame to I-frame and vice versa
%
arg    = tau+THETA;
carg   = cos(arg);
sarg   = sin(arg);
cdelta = cos(delta);
sdelta = sin(delta);
%
Civ(1,1) = -sdelta*carg;
Civ(2,1) = -sdelta*sarg;
Civ(3,1) =  cdelta;
Civ(1,2) = -sarg;
Civ(2,2) =  carg;
Civ(3,2) =  0;
Civ(1,3) = -cdelta*carg;
Civ(2,3) = -cdelta*sarg;
Civ(3,3) = -sdelta;
%
Cvi      = Civ';
%
% Transformtion I-frame to R-frame and vice versa
%
cangle = cos(THETA);
sangle = sin(THETA);
%
Cri      = zeros(3,3);
Cri(1,1) =  cangle;
Cri(1,2) =  sangle;
Cri(2,1) = -sangle;
Cri(2,2) =  cangle;
Cri(3,3) =  1;
%
Cir      = Cri';
%
% Transformation R-frame to V-frame and vice versa
%
Cvr = Cvi*Cir;
Crv = Cvr';
%
% Spherical (relative) to Cartesian (relative) position
%
cdelta = cos(delta);
Xr(1)  = r*cdelta*cos(tau);
Xr(2)  = r*cdelta*sin(tau);
Xr(3)  = r*sin(delta);
%
% Inertial Cartesian position
%
Xi = Cir*Xr';
%
% Spherical velocity to Cartesian V-frame (relative)
%
Vv(1) =  Vrel*cos(gamma)*cos(chi);
Vv(2) =  Vrel*cos(gamma)*sin(chi);
Vv(3) = -Vrel*sin(gamma);
Crv   = Cri*Civ;
Vr    = Crv*Vv';
Vi    = Cir*Vr + cross([0 0 omcb]',Xi);
%%
% INVERSE
%
Xr = Cri*Xi;
%
r     = norm(Xr);
h     = r - Re;
delta = asin(Xr(3)/r);
tau   = atan2(Xr(2),Xr(1));
%
% Cartesian velocity from I-frame to R-frame
%
Vxri = Vi(1) + omcb*Xi(2);
Vyri = Vi(2) - omcb*Xi(1);
%
Vr(1) = Cri(1,1)*Vxri + Cri(1,2)*Vyri;
Vr(2) = Cri(2,1)*Vxri + Cri(2,2)*Vyri;
Vr(3) = Vi(3);
%
% Cartesian velocity in V-frame (NED)
%
Vv    = Cvr*Vr;
V     =  norm(Vv);
gamma = -asin(Vv(3)/V);
chi   =  atan2(Vv(2),Vv(1));
