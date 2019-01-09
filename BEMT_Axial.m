#!/usr/bin/octave

clear; clc; clf;

% With collective ONLY
% No twist, coning
% Hover flight

% Specified Rotor and Aerodynamic Parameters
R = 8.0;
AR = 6.0;                 % Aspect ratio
c = R/AR;;
Nb = 2;
Om_rpm = 1250;            % in rpm
a = 5.74;                 % d_Cl/d_alpha in radians
rho = 1.2;                % in kg/m3
root_cut=0.05;
tip_cut=0.90;             % For accounting tip loss
vel_climb = 0.0;
coning=0.0;               % Not accounted for currently

theta=8*pi/180;

% Calculated Rotor Parameters
sol = Nb*c/(pi*R);
Om = Om_rpm*pi/30;
lam_c = vel_climb/(R*Om);

% Solver parameters 
nx = 50;                    % No. of stations along blade

% Initialization
% Equi-spaced
r_bar = linspace(root_cut,tip_cut,nx);
dr_bar = r_bar(2)-r_bar(1);

% Cosine spaced
%theta_cosine=linspace(0,pi,nx);
%r_bar = R*0.5-R*0.5*cos(theta_cosine);
%dr_bar=r_bar(1:end-1)-r_bar(2:end);
%dr_bar=-1*[dr_bar 0];

% atan spaced
%x_spacing=linspace(tan(-1),tan(1),nx);
%r_bar=R*0.5+R*0.5*atan(x_spacing);
%dr_bar=r_bar(2:end)-r_bar(1:end-1);
%dr_bar=[dr_bar 0];

const1 = sol*a/16; 
const_c = lam_c*0.5;

% for theta=degtorad(10:0.1:12)

% Inflow ratio from BEMT
lam = -(const1-const_c) + sqrt((const1-const_c)^2+2*const1*theta*r_bar);
alf=theta-lam./r_bar;

plot(r_bar,lam);
grid on;
xlabel('Normalized Radius');
ylabel('Inflow Ratio');

% Using Momentum theory
ct_vec = 4*lam.*(lam-lam_c).*r_bar.*dr_bar;
format long;
CT_MT = sum(ct_vec)

% Using BEMT
ct_vec = 0.5*sol*a.*dr_bar.*(r_bar.^2).*alf;
CT_BEMT = sum(ct_vec);
% plot(r_bar,ct_vec.*(rho*pi*R*R*(R*Om)^2))

Thrust = CT_MT*(rho*pi*R*R*(R*Om)^2)

% Sectional lift distribution
cl_vec = 2.0*ct_vec/(c/R*dr_bar);

% plot(radtodeg(theta),CT,'.')
% hold on
% end
% grid on
lam_i=lam;
% fid=fopen('Vz_OGE_BEMT.dat','w');
% fprintf(fid,'Variables = "r/R" "vz"\n');
% fprintf(fid,'Zone T= "Vz_OGE_BEMT"\n');
% dlmwrite('Vz_OGE_BEMT.dat',[r_bar'*R lam_i'*R*Om],'-append','delimiter',' ');
