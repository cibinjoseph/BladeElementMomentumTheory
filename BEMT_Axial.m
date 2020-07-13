clear; clc; clf;

% With collective ONLY
% Axial flight condition

% Rotor parameters
R = 0.3048;                % in m
Nb = 2;                    % blade number
Om_rpm = 6000;             % in rpm
vel_climb = 18.0;          % in m/s

CLa = 4.88; %2*pi;         % d_Cl/d_alpha in radians
CL0 = 0.470;

% Selective parameters
% *** COMMENT OUT UNUSED VARIABLES ***
% -- AR or chord --
% AR = 6;                % Aspect ratio
% c = 0;
cByR_dist = [  % as a two-col matrix [r/R  c/R]
  0.25000   0.08
  0.29137   0.09842
  0.32614   0.1119
  0.35372   0.12106
  0.38129   0.1295
  0.41007   0.13746
  0.43645   0.143
  0.46763   0.14712
  0.50360   0.15028
  0.55276   0.15176
  0.60911   0.15062
  0.69305   0.1459
  0.76259   0.13972
  0.82734   0.13208
  0.90048   0.12086
  1.00000   0.10054
  ];
% -- theta distribution or constant angle --
% theta_deg = 8;
theta_deg_dist = [  % as a two-col matrix [r/R theta]
  0.25	32.56026
  0.28037	30.10193
  0.3007	28.30318
  0.33182	25.96551
  0.36413	23.8082
  0.39287	22.1909
  0.42521	20.51407
  0.46234	18.95808
  0.49709	17.7621
  0.53783	16.4469
  0.59416	15.01391
  0.64569	13.82045
  0.70084	12.86777
  0.75599	12.03521
  0.82912	11.20535
  0.89987	10.49525
  1	10.21008
  ]; 
% ------

% Environment parameters
rho = 1.2;                % in kg/m3

% Solver parameters 
nx = 50;                  % No. of stations along blade

% Feature switches
spacing_switch = 1;       % [1]Equispaced [2]Cosine [3]aTan
prandtlTipLoss_switch = 1;

% For accounting tip loss
% (change accordingly when prandtlTipLoss_switch is 0)
root_cut = 0.25;            % r/R
tip_cut = 0.999;             % r/R

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

% Calculated Rotor Parameters
if (exist('AR') == 1)
  c = ones(1, nx) * R/AR;
elseif (exist('c') == 1) 
  c = ones(1, nx) * c;
elseif (exist('cByR_dist') == 1)
  cByR_poly = polyfit(cByR_dist(:, 1), cByR_dist(:, 2), 4);
  cByR = polyval(cByR_poly, r_bar);
  c = cByR*R;
end

if (exist('theta_deg') == 1)
  theta_deg = theta_deg*ones(1,nx);
elseif (exist('theta_deg_dist') == 1)
  theta_poly = polyfit(theta_deg_dist(:, 1), theta_deg_dist(:, 2), 4);
  theta_deg = polyval(theta_poly, r_bar);
  plot(theta_deg_dist(:,1), theta_deg_dist(:,2), 'bo'); hold on;
  plot(r_bar, theta_deg,'r-'); hold off;
end

alphaCamber = CL0/CLa;
theta = theta_deg*pi/180 + alphaCamber;  % Corrected for cambered airfoils
sol = Nb*c/(pi*R);
Om = Om_rpm*pi/30;        % in rad per sec
lam_climb = vel_climb/(R*Om);

% Inflow computation
switch prandtlTipLoss_switch

  case 0
    % Inflow ratio from BEMT
    const1 = sol*CLa/16; 
    const_climb = lam_climb*0.5;
    lam = -(const1-const_climb) ...
      + sqrt((const1-const_climb).^2+2*const1.*theta.*r_bar);
    phi=lam./r_bar;
    alf=theta-phi;
    prandtl_F = ones(size(r_bar));

  case 1
    prandtl_F = ones(size(r_bar));  % Initial value
    prandtl_residual = 1.0;
    const1 = sol*CLa/16; 
    const_climb = lam_climb*0.5;
    counter = 1;

    while ( (prandtl_residual>0.001) && (counter<21) )
      lam = -(const1./prandtl_F-const_climb) ...
        + sqrt((const1./prandtl_F-const_climb).^2 ...
        + 2*const1./prandtl_F.*theta.*r_bar);
      prandtl_Fsmall = 0.5*Nb*(1.0-r_bar)./lam;

      prandtl_F_prev = prandtl_F;
      prandtl_F = (2/pi)*acos(exp(-prandtl_Fsmall));
      prandtl_residual = abs(norm(prandtl_F-prandtl_F_prev));
      counter = counter+1;
    end

    if (counter==21)
      warning('Prandtl tip-loss factor failed to converge');
    end

    phi=lam./r_bar;
    alf=theta-phi;

end

% Using Momentum theory
ct_vec = prandtl_F.*4.*lam.*(lam-lam_climb).*r_bar.*dr_bar;
format long;
CT_MT = sum(ct_vec);

% Using BEMT
ct_vec = 0.5*sol*CLa.*dr_bar.*(r_bar.^2).*alf;
CT_BEMT = sum(ct_vec);

Thrust = CT_MT*(rho*pi*R*R*(R*Om)^2);

% Check between MT and BEMT
if (CT_MT-CT_BEMT)>eps
  warning('Warning: Discrepancy between CT calculated using BEMT and MT')
end

% Sectional lift distribution
cl_vec = 2.0*ct_vec/(c./R.*dr_bar);

% Results
% fprintf('\nColl. pitch (deg) = %d\n',theta_deg);
% fprintf('Solidity = %d\n\n',sol);
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
plot(r_bar,(alf-alphaCamber)*180/pi,'k');
grid on;
xlabel('r/R');
ylabel('Alpha (deg)');

subplot(2,2,4);
plot(r_bar,alf*CLa,'k');
grid on;
xlabel('r/R');
ylabel('Sectional CL');

% Write to file
% fid=fopen('Vz_OGE_BEMT.dat','w');
% fprintf(fid,'Variables = "r/R" "vz"\n');
% fprintf(fid,'Zone T= "Vz_OGE_BEMT"\n');
% dlmwrite('Vz_OGE_BEMT.dat',[r_bar'*R lam'*R*Om],'-append','delimiter',' ');
