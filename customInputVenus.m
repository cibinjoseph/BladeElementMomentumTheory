% With collective ONLY
% Axial flight condition

% Rotor parameters
R = 1.0;                % in m
Nb = 4;                    % blade number
Om_rpm = 74.18;             % in rpm
vel_climb = 0.0;          % in m/s

CLa = 2*pi;         % d_Cl/d_alpha in radians
alpha0_deg = 0.0;

% Selective parameters
% *** COMMENT OUT UNUSED VARIABLES ***
% -- AR or chord --
% AR = 6;                % Aspect ratio
c = 0.1;
% cByR_dist = [  % as a two-col matrix [r/R  c/R]
%   0.2	0.116
%   0.3	0.109
%   0.4	0.101
%   0.5	0.094
%   0.6	0.087
%   0.7	0.080
%   0.8	0.072
%   0.9	0.065
%   1.0	0.058
%   ];
% -- theta distribution or constant angle --
% theta_deg = 11.348;
theta_deg_dist = [  % as a two-col matrix [r/R theta]
  0.3	11.348
  1.0	11.348-3.5362
  ]; 
% ------

% Environment parameters
rho = 0.2054;                % in kg/m3

% Solver parameters 
nx = 50;                  % No. of stations along blade

% Feature switches
spacing_switch = 1;       % [1]Equispaced [2]Cosine [3]aTan
prandtlTipLoss_switch = 1;

% For accounting tip loss
% (change accordingly when prandtlTipLoss_switch is 0)
root_cut = 0.3;            % r/R
tip_cut = 0.999;             % r/R
