%% Simulation Data
simu = simulationClass();               % Initialize Simulation Class
simu.simMechanicsFile = 'OSWEC.slx';    % Specify Simulink Model File
simu.mode = 'normal';                   % Specify Simulation Mode ('normal','accelerator','rapid-accelerator')
simu.explorer='off';                     % Turn SimMechanics Explorer (on/off)
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = 50;                     % Wave Ramp Time [s]
simu.endTime=200;                       % Simulation End Time [s]        
simu.solver = 'ode45';                   % simu.solver = 'ode4' for fixed step & simu.solver = 'ode45' for variable step 
simu.dt = 0.01;                          % Simulation Time-Step [s]
simu.CITime = 30;                       % Specify CI Time [s]

%% Wave Information
% % noWaveCIC, no waves with radiation CIC  
% waves = waveClass('noWaveCIC');       % Initialize Wave Class and Specify Type  

% % Regular Waves 
% waves = waveClass('regular');           % Initialize Wave Class and Specify Type                                 
% waves.H = 2.5;                          % Wave Height [m]
% waves.T = 8;                            % Wave Period [s]

% Irregular Waves using PM Spectrum with Directionality 
waves = waveClass('irregular');         % Initialize Wave Class and Specify Type
waves.H = 2.5;                          % Significant Wave Height [m]
waves.T = 8;                            % Peak Period [s]
waves.spectrumType = 'PM';              % Specify Spectrum Type
waves.phaseSeed = 8;
%waves.waveDir = [0,30,90];              % Wave Directionality [deg]
%waves.waveSpread = [0.1,0.2,0.7];       % Wave Directional Spreading [%}

% % Irregular Waves with imported spectrum
% waves = waveClass('spectrumImport');        % Create the Wave Variable and Specify Type
% waves.spectrumDataFile = 'spectrumData.mat';  %Name of User-Defined Spectrum File [:,2] = [f, Sf]

% % Waves with imported wave elevation time-history  
% waves = waveClass('etaImport');         % Create the Wave Variable and Specify Type
% waves.etaDataFile = 'etaData.mat'; % Name of User-Defined Time-Series File [:,2] = [time, eta]

ptosim = ptoSimClass('Non-Compressible Fluid Hydraulic');

% Rotary to Linear Adjustable Rod
% same sizing as PTO-Sim OSWEC example
ptosim.motionMechanism.crank = 3;
ptosim.motionMechanism.offset = 1.3;
ptosim.motionMechanism.rodInit = 5;

%% Body Data
% Flap
body(1) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Flap
body(1).geometryFile = 'geometry/flap.stl';     % Geometry File
body(1).mass = 127000;                          % User-Defined mass [kg]
body(1).momOfInertia = [1.85e6 1.85e6 1.85e6];  % Moment of Inertia [kg-m^2]

% Base
body(2) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Base
body(2).geometryFile = 'geometry/base.stl';     % Geometry File
body(2).mass = 'fixed';                         % Creates Fixed Body

%% PTO and Constraint Parameters
% Fixed
constraint(1)= constraintClass('Constraint1');  % Initialize ConstraintClass for Constraint1
constraint(1).loc = [0 0 -10];                  % Constraint Location [m]

%kp = 2.9897e7;
%ki = 3.2542e7;

% Rotational PTO
pto(1) = ptoClass('PTO1');                      % Initialize ptoClass for PTO1
%pto(1).k = kp;                                 % PTO Stiffness Coeff [Nm/rad]
%pto(1).c = ki;                               % PTO Damping Coeff [Nsm/rad]
pto(1).loc = [0 0 -8.9];                        % PTO Location [m]

% grid seach - work in - regular waves
%pto(1).c = 1.9e07;                            % PTO Damping Coeff [Nsm/rad]
%pto(1).k = 3.3e7;                                   % PTO Stiffness Coeff [Nm/rad]

% grid seach - work in irregular waves
%pto(1).c = 3.1e07;                             % PTO Damping Coeff [Nsm/rad]
%pto(1).k = 3.1e7;                                   % PTO Stiffness Coeff [Nm/rad]


% grid seach - EHA work out
%pto(1).c = 2.875e07;     % pto(1).c = 3.1e07;                             % PTO Damping Coeff [Nsm/rad]
%pto(1).k = 2.875e7;       %  pto(1).k = 2.9e7;                                   % PTO Stiffness Coeff [Nm/rad]

% grid seach - EHA work out - irregular waves
%pto(1).c = 1.6875e7; % Scale = 130 
pto(1).c = 2.0833e7;       %  pto(1).c = 4.1e07;                             % PTO Damping Coeff [Nsm/rad]
%pto(1).k = 3.2708e7; % Scale = 130 
pto(1).k = 3.6667e7;          %pto(1).k = 2.1e7;                                   % PTO Stiffness Coeff [Nm/rad]


%Perry's - from maximum power of linear system
% pto(1).c = 2.0188e7;                             % PTO Damping Coeff [Nsm/rad]
% pto(1).k = 3.8278e7;                                 % PTO Stiffness Coeff [Nm/rad]



%% Rotary to Linear Adjustable Rod
%ptosim.motionMechanism.crank = 3;
%ptosim.motionMechanism.offset = 1.3;
%ptosim.motionMechanism.rodInit = 5;
