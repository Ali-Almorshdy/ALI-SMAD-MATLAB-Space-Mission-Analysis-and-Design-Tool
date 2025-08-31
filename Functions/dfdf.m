function [dydt]=dfdf(t,y,JD0)
x1 = y(1);
y1 = y(2);
z1 = y(3);
ae = 6378;             
GM = 398600.44189;
% GMS=132.712e9;
r1 = norm(y(1:3));
wE = [ 0 0 7.2921159e-5]';
A =20;
alt = r1 - ae;
rho = Aliatmosphere(alt);
V = y(4:6)' ;
Vrel = V - cross(wE,y(1:3)');
vrel = norm(Vrel); %Speed relative to the atmosphere (km/s)
uv = Vrel/vrel; %Relative velocity unit vector
ap = -2.2*A/500*rho*... %Acceleration due to drag (m/s^2)
    (1000*vrel)^2/2*uv; %(converting units of vrel from km/s to m/s)
JD = JD0 + t/86400;
[~, ~, ~,r_S] = solar_position(JD);
SR=r_S-y(1:3);
SRN=norm(SR);
usr=SR/SRN;
%  rsubdot = -GMS*SR/SRN^3;
nu = los(y(1:3), r_S);
c = 2.998e8; %Speed of light (m/s)
S = 1367; %Solar constant (W/m^2)
%...Satellite data:
CR = 1; %Radiation pressure codfficient
m = 500; %Mass (kg)
As = 20; %Frontal area (m^2);
pSR = (nu*(S/c)*CR*As/m/1000)*usr;
J2 =  1.0826353865466185e-03;
ax = -GM/r1^3*x1*(1 ...
   -J2*(3/2)*(ae/r1)^2*(5*(z1/r1)^2-1));
ay = y1/x1*ax;
az = -GM/r1^3*z1*(1 ...
   +J2*(3/2)*(ae/r1)^2*(3-5*(z1/r1)^2));
dydt = [y(4:6);[ax;ay;az]+ap'/1000+pSR];




end