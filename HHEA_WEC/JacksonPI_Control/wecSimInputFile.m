%% Simulation Data
simu = simulationClass();               % Initialize simulationClass
simu.simMechanicsFile = 'OSWEC.slx';        % Simulink Model File
simu.startTime = 0;                     % Simulation Start Time [s]
simu.rampTime = 100;                    % Wave Ramp Time [s]
simu.endTime=400;                       % Simulation End Time [s]        
simu.dt = 0.1;                          % Simulation Time-Step [s]

%% Wave Information
% Regular Waves  
waves = waveClass('regular');              % Initialize waveClass                                
waves.T = 2.5;                          % Wave Period [s]
waves.H = 8;                          % Wave Height [m]

%% Body Data
% Flap
body(1) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Flap
body(1).geometryFile = 'geometry/flap.stl';     % Geometry File
body(1).mass = 127000;                          % Mass [kg]
body(1).momOfInertia = [1.85e6 1.85e6 1.85e6];  % Moment of Inertia [kg-m^2]
%body(1).initDisp = struct( 'initLinDisp', [0 0 0], 'initAngularDispAxis', [0 1 0], 'initAngularDispAngle', -.13);  % Initial condition- where is the body's center of mass at t = 0 [m] [x y z]

% Base
body(2) = bodyClass('hydroData/oswec.h5');      % Initialize bodyClass for Base
body(2).geometryFile = 'geometry/base.stl';     % Location of Geomtry File
body(2).mass = 'fixed';                         % Creates Fixed Body

%% PTO and Constraint Parameters
% Fixed
constraint(1)= constraintClass('Constraint1');  % Initialize constraintClass for Constraint1
constraint(1).loc = [0 0 -10];                  % Constraint Location [m]

% Rotational PTO
pto(1) = ptoClass('PTO1');                      % Initialize ptoClass for PTO1
pto(1).k = 0;                                	% PTO Stiffness Coeff [Nm/rad]
pto(1).c = 10;                                   % PTO Damping Coeff [Nsm/rad]
pto(1).loc = [0 0 -8.9];                        % PTO Location [m]