clear; clc; clf;

% With collective ONLY
% Axial flight condition

isCustomInput = true;

if isCustomInput
  fprintf('CUSTOM input files used\n')
  customInput;
else
  fprintf('DEFAULT input files used\n')

  % Rotor parameters
  R = 0.518;                % in m
  Nb = 2;                    % blade number
  Om_rpm = 3000;             % in rpm
  vel_climb = 2.0;          % in m/s

  CLa = 2*pi;         % d_Cl/d_alpha in radians
  alpha0_deg = 0.0;

  % Selective parameters
  % *** COMMENT OUT UNUSED VARIABLES ***
  % -- AR or chord --
  % AR = 6;                % Aspect ratio
  % c = 0;
  cByR_dist = [  % as a two-col matrix [r/R  c/R]
    0.2	0.116
    0.3	0.109
    0.4	0.101
    0.5	0.094
    0.6	0.087
    0.7	0.080
    0.8	0.072
    0.9	0.065
    1.0	0.058
    ];
  % -- theta distribution or constant angle --
  % theta_deg = 8;
  theta_deg_dist = [  % as a two-col matrix [r/R theta]
    0.2	14.97
    0.3	13.97
    0.4	12.97
    0.5	11.97
    0.6	10.97
    0.7	9.97
    0.8	8.97
    0.9	7.97
    1.0	6.97
    ]; 
  % ------

  % Environment parameters
  rho = 0.022;                % in kg/m3

  % Solver parameters 
  nx = 50;                  % No. of stations along blade

  % Feature switches
  spacing_switch = 1;       % [1]Equispaced [2]Cosine [3]aTan
  prandtlTipLoss_switch = 1;

  % For accounting tip loss
  % (change accordingly when prandtlTipLoss_switch is 0)
  root_cut = 0.20;            % r/R
  tip_cut = 0.999;             % r/R
end

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
  if (size(theta_deg_dist, 1) == 2)
    theta_poly = polyfit(theta_deg_dist(:, 1), theta_deg_dist(:, 2), 1);
  else
    theta_poly = polyfit(theta_deg_dist(:, 1), theta_deg_dist(:, 2), 4);
  end
  theta_deg = polyval(theta_poly, r_bar);
  plot(theta_deg_dist(:,1), theta_deg_dist(:,2), 'bo'); hold on;
  plot(r_bar, theta_deg,'r-'); hold off;
end
% disp('DBG')
% disp(polyval(theta_poly, 0.75))
% disp('DBG END')
% keyboard

alphaCamber = -1.0*alpha0_deg*pi/180.0;
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
lamMean = 2.0*trapz(r_bar, lam.*r_bar)/(1-root_cut*root_cut);
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
cl_vec = 2.0*ct_vec./(sol.*r_bar.*r_bar.*dr_bar);
[clMax, imx] = max(cl_vec);
gammaMax = 0.5*clMax*c(imx)*(R*dr_bar(imx))*norm([r_bar(imx) lam(imx)])*R*Om;

% Results
% fprintf('\nColl. pitch (deg) = %d\n',theta_deg);
% fprintf('Solidity = %d\n\n',sol);
fprintf('CT           = %d\n', CT_BEMT);
fprintf('inflow ratio = %d\n', lamMean);
fprintf('Thrust (N)   = %d\n', Thrust);
fprintf('inflow (m/s) = %d\n', lamMean*R*Om);
fprintf('CL max       = %d\n', clMax);
fprintf('Gamma max    = %d\n', gammaMax);

% Generate plots
subplot(2,2,1);
plot(r_bar,lam,'k');
grid on;
xlabel('r/R');
ylabel('Inflow Ratio');

subplot(2,2,2);
plot(r_bar,ct_vec/Nb,'k');
grid on;
xlabel('r/R');
ylabel('Sectional CT');

subplot(2,2,3);
plot(r_bar,(alf-alphaCamber)*180/pi,'k');
grid on;
xlabel('r/R');
ylabel('Alpha (deg)');

subplot(2,2,4);
plot(r_bar,(cl_vec),'k');
grid on;
xlabel('r/R');
ylabel('Sectional CL');

% Write to file
% fid=fopen('Vz_OGE_BEMT.dat','w');
% fprintf(fid,'Variables = "r/R" "vz"\n');
% fprintf(fid,'Zone T= "Vz_OGE_BEMT"\n');
% dlmwrite('Vz_OGE_BEMT.dat',[r_bar'*R lam'*R*Om],'-append','delimiter',' ');
