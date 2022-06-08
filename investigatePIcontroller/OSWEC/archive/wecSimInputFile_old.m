%% Simulation Data
simu = simulationClass();              
simu.simMechanicsFile = 'OSWEC_basic_PI.slx'; % Specify Simulink Model File with PTO-Sim
simu.mode='normal';
simu.explorer = 'off';

simu.rampTime = 50;
simu.endTime = 200;  
simu.dt = 0.01;
simu.dtOut = 10*simu.dt;            % Variable output time step. Good for Paraview speed
simu.CITime = 30;                 % Convolution Integral time. (default=60, using 30 should be faster)
simu.solver = 'ode4';
% simu.ssCalc = 1;                    % Do state space calculation instead of CIC, can't do with irregular?
simu.rho = 1025;

% % Body to Body effects
simu.b2b = 0;                       % Turn B2B interactions 'on'

% % Paraview options (turn off for MCR)
% %    See files in \vtk for Paraview writing progress
% simu.paraview = 1;             	% Saves data to *.vtp for Paraview
% simu.pressureDis = 1;             % Saves pressures data for Paraview. Only for nonlinear hydro

% % Nonlinear Hydro effects
simu.nlHydro = 2;
simu.dtNL = 1*simu.dt;

%% Wave Information
%Irregular Waves using PM Spectrum
waves = waveClass('irregular');
waves.H = 2.5;
waves.T = 8;
waves.spectrumType = 'BS';
waves.phaseSeed=1;

%% Body Data
% Flap
body(1) = bodyClass('../hydroData/oswec.h5');   
body(1).geometryFile = '../geometry/flap.stl';    
body(1).mass = 127000;
body(1).momOfInertia = [1.85e6 1.85e6 1.85e6];
body(1).linearDamping = [0, 0, 0, 0, 1*10^7, 0];    % Specify damping on body

% Base
body(2) = bodyClass('../hydroData/oswec.h5');   
body(2).geometryFile = '../geometry/base.stl';    
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

%\% Mooring
% Moordyn
% mooring(1) = mooringClass('mooring');       % Initialize mooringClass
% mooring(1).moorDynLines = 6;                % Specify number of lines
% mooring(1).moorDynNodes(1:3) = 16;          % Specify number of nodes per line
% mooring(1).moorDynNodes(4:6) = 6;           % Specify number of nodes per line
% mooring(1).initDisp.initLinDisp = [0 0 -0.21];      % Initial Displacement

% Mooring Matrix
% mooring(1) = mooringClass('mooring');       % Initialize mooringClass
% mooring(1).matrix.k = zeros(6,6);
% mooring(1).matrix.k(1,1) = 1e5;
% mooring(1).matrix.c = zeros(6,6);
% mooring(1).matrix.preTension = zeros(1,6);

%% Controller Parameters
controller = struct();
controller.name = 'PI';
controller.Kp = 0;
controller.Ki = 0;
