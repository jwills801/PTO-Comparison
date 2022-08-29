% See So et al. 2015 PTO-Sim paper: 
% doi: 10.1109/PESGM.2015.7285735
% https://ieeexplore.ieee.org/abstract/document/7285735

%% Simulation Data
simu = simulationClass();              
%simu.simMechanicsFile = 'OSWEC_Hydraulic_PTO.slx'; % Specify Simulink Model File with PTO-Sim
simu.simMechanicsFile = 'OSWEC_Constant_DeltaP.slx'; % Specify Simulink Model File
simu.mode = 'normal';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
%simu.mode = 'rapid-accelerator';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer = 'off';                   % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = 50;                    % Wave Ramp Time [s]
simu.endTime = 2000;                     % Simulation End Time [s]        
simu.solver = 'ode4';                   % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.01;                         % Simulation Time-Step [s]
%simu.rho = 1025;

simu.CITime = 30;                       % Specify CI Time for OSWEC, max 30s
%simu.b2b = 0;                           % doesn't matter bc base is fixed
% simu.ssCalc = 1;                        % Turn state space on

% simu.nlHydro = 2;                       % giving errors with irregular wave case, dont use for now 
% simu.dtNL = 1*simu.dt;
% --> irregular+nlHydro having issues?? check this with other t_osu issue case as well

%% Wave Information
% % Regular wave case (A)
% waves = waveClass('regular');           % Initialize Wave Class and Specify Type
% waves.H = 2.5;                          % Significant Wave Height [m]
% waves.T = 8;                            % Peak Period [s]

% Irregular wave case (B) using PM Spectrum. No directionality
waves = waveClass('irregular');         % Initialize Wave Class and Specify Type
waves.H = 2.5;                          % Significant Wave Height [m]
waves.T = 8;                            % Peak Period [s]
waves.spectrumType = 'PM';              % PM Spectrum Type
waves.phaseSeed = 8;

%% Body Data
% Flap
body(1) = bodyClass('../hydroData/oswec.h5');      % Initialize bodyClass for Flap
body(1).geometryFile = '../geometry/flap.stl';     % Geometry File
body(1).mass = 127000;                          % User-Defined mass [kg]
body(1).momOfInertia = [1.85e6 1.85e6 1.85e6];  % Moment of Inertia [kg-m^2]
body(1).linearDamping = zeros(6);               % Specify damping on body
% body(1).linearDamping(5,5) = 1e7; % set pitch linear damping for Kp
body(1).viscDrag.cd = [1 1 1 1 1 1]*1.8; % drag coefficient

% oswec dims: 
x=1.8; y=9; z=10.6;

% giorgio paper (https://doi.org/10.1016/j.arcontrol.2015.09.006):
%     dimensions: x=1, y=30, z=15, rho_b=250
%     drag coefficient: cd=1.9,1.8
%     drag moment_y = -B_v2 * w * |w| (w=dx/dt(5)) 

body(1).viscDrag.characteristicArea = [y*z x*z x*y 0 (y*z^4)/4 0];
% in WEC-Sim: F_viscDrag = 0.5*rho.*cd.*ca * v * |v|
% in giorgio paper F_viscDrag = -Bv2 * v * |v|, Bv2 = 1/8 rho cd w h^4
%  1/2 rho cd (ca)        v |v|
%  1/2 rho cd (1/4 w h^4) v |v|
%      Solving for ca in pitch gives: (w*h^4)/4
% ca_4 and ca_6 dont matter bc its not moving in roll/yaw

clear x y z
% Base
body(2) = bodyClass('../hydroData/oswec.h5');   
body(2).geometryFile = '../geometry/base_shift.stl';    
body(2).mass = 'fixed';   

%% PTO and Constraint Parameters
% Fixed Constraint
constraint(1)= constraintClass('Constraint1');  
constraint(1).loc = [0 0 -10];   

% Rotational PTO
pto(1) = ptoClass('PTO1');                      % Initialize ptoClass for PTO1
pto(1).k = 0;                                   % PTO Stiffness Coeff [Nm/rad]
pto(1).c = 0;                                   % PTO Damping Coeff [Nsm/rad]
pto(1).loc = [0 0 -8.9];                        % PTO Location [m]


