%% Non-Compressible Fluid Hydraulic PTO-Sim
% See So et al. 2015 PTO-Sim paper: 
% doi: 10.1109/PESGM.2015.7285735
% https://ieeexplore.ieee.org/abstract/document/7285735

ptosim = ptoSimClass('Non-Compressible Fluid Hydraulic');

%% Piston 

ptosim.pistonNCF.topA = 0.2382;    % Top piston area [m^2]
ptosim.pistonNCF.botA = 0.2382;    % Bottom piston area [m^2]

%ptosim.pistonNCF.topA = 0.2;    % Top piston area [m^2]
%ptosim.pistonNCF.botA = 0.2;    % Bottom piston area [m^2]


%% Low Pressure Accumulator
% ptosim.accumulator(2).VI0 = 6;                                                           % Initial volume [m^3]
% ptosim.accumulator(2).pIrated = 16e6;                                                    % Rated working pressure
% ptosim.accumulator(2).pIupper_limit = (4/3)*ptosim.accumulator(2).pIrated;               % Upper working pressure
% ptosim.accumulator(2).pIlower_limit = (0.5)*ptosim.accumulator(2).pIupper_limit;         % Lower working pressure
% ptosim.accumulator(2).pIprecharge = 0.9*ptosim.accumulator(2).pIlower_limit;             % Precharge pressure
% ptosim.accumulator(2).VImax = ptosim.accumulator(2).VI0*(1-(ptosim.accumulator(2).pIprecharge/ptosim.accumulator(2).pIupper_limit)^(1/1.4));
% ptosim.accumulator(2).VImin = ptosim.accumulator(2).VI0*(1-(ptosim.accumulator(2).pIprecharge/ptosim.accumulator(2).pIlower_limit)^(1/1.4));
% ptosim.accumulator(2).VIeq = ptosim.accumulator(2).VImax;
% ptosim.accumulator(2).pIeq = ptosim.accumulator(2).pIprecharge/(1-ptosim.accumulator(2).VIeq/ptosim.accumulator(2).VI0)^(1.4);

%ptosim.accumulator(2).VI0 = X(1)/2; % Initial volume [m^3] - assuming that the initial volume is half full of fluid
%ptosim.accumulator(2).pIprecharge = X(3) *2^1.2; % Initial pressure [Pa] - assuming that the initial volume is half full of fluid
ptosim.accumulator(2).del_p_r = NaN;                                         
ptosim.accumulator(2).pIrated = NaN;
ptosim.accumulator(2).pIeq = NaN;
ptosim.accumulator(2).pIlower_limit = NaN;
ptosim.accumulator(2).pIupper_limit = NaN;
ptosim.accumulator(2).VIeq = 0;
ptosim.accumulator(2).VImax = NaN;
ptosim.accumulator(2).VImin = NaN;

%% High Pressure Accumulator
% ptosim.accumulator(1).VI0 = 8.5;                                                                   % Initial volume                          
% ptosim.accumulator(1).del_p_r = 10e6;                                         
% ptosim.accumulator(1).pIrated = ptosim.accumulator(1).del_p_r + ptosim.accumulator(2).pIrated;     % Rated working pressure
% ptosim.accumulator(1).pIeq = ptosim.accumulator(2).pIeq + 1e7;
% ptosim.accumulator(1).pIlower_limit = ptosim.accumulator(1).pIeq;
% ptosim.accumulator(1).pIupper_limit = 1.5*ptosim.accumulator(1).pIlower_limit;
% ptosim.accumulator(1).pIprecharge = 0.4*ptosim.accumulator(1).pIlower_limit;
% ptosim.accumulator(1).VIeq = ptosim.accumulator(1).VI0*(1-(ptosim.accumulator(1).pIprecharge/ptosim.accumulator(1).pIeq)^(1/1.4));
% ptosim.accumulator(1).VImax = ptosim.accumulator(1).VI0*(1-(ptosim.accumulator(1).pIprecharge/ptosim.accumulator(1).pIupper_limit)^(1/1.4));
% ptosim.accumulator(1).VImin = ptosim.accumulator(1).VI0*(1-(ptosim.accumulator(1).pIprecharge/ptosim.accumulator(1).pIlower_limit)^(1/1.4));

%ptosim.accumulator(1).VI0 = X(2)/2; % Initial volume [m^3] - assuming that the initial volume is half full of fluid
%ptosim.accumulator(1).pIprecharge = X(4) *2^1.2; % Initial pressure [Pa] - assuming that the initial volume is half full of fluid
ptosim.accumulator(1).del_p_r = NaN;                                         
ptosim.accumulator(1).pIrated = NaN;
ptosim.accumulator(1).pIeq = NaN;
ptosim.accumulator(1).pIlower_limit = NaN;
ptosim.accumulator(1).pIupper_limit = NaN;
ptosim.accumulator(1).VIeq = 0;
ptosim.accumulator(1).VImax = NaN;
ptosim.accumulator(1).VImin = NaN;

%% Hydraulic Motor
ptosim.hydraulicMotor.angVelInit = 0;                       % Initial speed
ptosim.hydraulicMotor.J = 20;                               % Total moment of inertia (motor & generator) [kg-m^2]
ptosim.hydraulicMotor.fric = 0.05;                          % Fricton [kg-m^2/s]


%% Lookup Table Generator
load('../../AP_107cc_Sept16.mat')
Pmap = P1_Mapping;
Qmap = Q1_Act_Mapping;
Tmap = T1_Act_Mapping;
wmap = w1rad_Mapping;
w = 3000*2*pi/60; % Angular Velocity [rad/s]
% if length(waves.type) == length('irregular')
%     d_scale = 250/107; % irregular waves
% else
%     d_scale = 1200/107; % regular waves
% end
%d_scale = X(5)/107; % how much bigger is this motor from the 107cc motor that the map was made for


%% Rotary to Linear Adjustable Rod
ptosim.motionMechanism.crank = 3;
ptosim.motionMechanism.offset = 1.3;
ptosim.motionMechanism.rodInit = 5;
