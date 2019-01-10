clear; clc; clf;

% With collective ONLY
% No twist, coning
% Hover flight condition

% Rotor parameters
R = 1.143;                % in m
AR = 6.0;                 % Aspect ratio
Nb = 2;
Om_rpm = 1250;            % in rpm
a = 5.74;                 % d_Cl/d_alpha in radians
vel_climb = 0.0;          % in m/s
theta_deg=8;

% For accounting tip loss
root_cut=0.05;            % r/R
tip_cut=0.90;             % r/R

% Environment parameters
rho = 1.2;                % in kg/m3

% Solver parameters 
nx = 50;                  % No. of stations along blade

% Feature switches
spacing_switch = 1;       % [1]Equispaced [2]Cosine [3]aTan
prandtlTipLoss_switch = 1;

% Calculated Rotor Parameters
c = R/AR;
sol = Nb*c/(pi*R);
theta=theta_deg*pi/180;
Om = Om_rpm*pi/30;        % in rad per sec
lam_climb = vel_climb/(R*Om);

% Spacing blade stations
switch spacing_switch

case 1
  % Equi-spaced
  r_bar = linspace(root_cut,tip_cut,nx);
  dr_bar = (r_bar(2)-r_bar(1))*ones(size(r_bar));

case 2
  % Cosine spaced
  theta_cosine=linspace(0,pi,nx);
  r_bar = R*0.5-R*0.5*cos(theta_cosine);
  dr_bar=r_bar(2:end)-r_bar(1:end-1);
  dr_bar=[dr_bar 0];

case 3
  % atan spaced
  x_spacing=linspace(tan(-1),tan(1),nx);
  r_bar=R*0.5+R*0.5*atan(x_spacing);
  dr_bar=r_bar(2:end)-r_bar(1:end-1);
  dr_bar=[dr_bar 0];

end

% Inflow computation
switch prandtlTipLoss_switch

case 0
  % Inflow ratio from BEMT
  const1 = sol*a/16; 
  const_c = lam_climb*0.5;
  lam = -(const1-const_c) + sqrt((const1-const_c)^2+2*const1*theta*r_bar);
  phi=lam./r_bar;
  alf=theta-phi;
  prandtl_F = 1;

case 1
  prandtl_F = 1;  % Initial value
  const1 = sol*a/16; 
  const_c = lam_climb*0.5;

  lam = -(const1./prandtl_F-const_c) + sqrt((const1./prandtl_F-const_c).^2+2*const1./prandtl_F*theta.*r_bar);
  phi=lam./r_bar;
  alf=theta-phi;

end

% Using Momentum theory
ct_vec = prandtl_F.*4*lam.*(lam-lam_climb).*r_bar.*dr_bar;
format long;
CT_MT = sum(ct_vec);

% Using BEMT
ct_vec = 0.5*sol*a.*dr_bar.*(r_bar.^2).*alf;
CT_BEMT = sum(ct_vec);

Thrust = CT_MT*(rho*pi*R*R*(R*Om)^2);

% Check between MT and BEMT
if (CT_MT-CT_BEMT)>eps
  warning('Warning: Discrepancy between CT calculated using BEMT and MT')
end

% Sectional lift distribution
cl_vec = 2.0*ct_vec/(c/R*dr_bar);

% Results
fprintf('\nColl. pitch (deg) = %d\n',theta_deg);
fprintf('Solidity = %d\n\n',sol);
fprintf('CT = %12.6f\n',CT_BEMT);
fprintf('Thrust (N) = %d\n',Thrust);

% Generate plots
subplot(2,2,1);
plot(r_bar,lam,'k');
grid on;
xlabel('r/R');
ylabel('Inflow Ratio');

subplot(2,2,2);
plot(r_bar,ct_vec,'k');
grid on;
xlabel('r/R');
ylabel('Sectional CT');

subplot(2,2,3);
plot(r_bar,alf*180/pi,'k');
grid on;
xlabel('r/R');
ylabel('Alpha (deg)');

subplot(2,2,4);
plot(r_bar,alf*a,'k');
grid on;
xlabel('r/R');
ylabel('Sectional CL');

% Write to file
%lam_i=lam;
% fid=fopen('Vz_OGE_BEMT.dat','w');
% fprintf(fid,'Variables = "r/R" "vz"\n');
% fprintf(fid,'Zone T= "Vz_OGE_BEMT"\n');
% dlmwrite('Vz_OGE_BEMT.dat',[r_bar'*R lam_i'*R*Om],'-append','delimiter',' ');
