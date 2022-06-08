% Script used to compare the power output of the PTO-Sim hydraulic PTO and 
% the HHEA system. 
% The HHEA uses the output of a PI controller used as a PTO in wec-sim. 
% Controller output contains the force and velocity of an actuator
% 
% See So et al. 2015 PTO-Sim paper: 
% doi: 10.1109/PESGM.2015.7285735
% https://ieeexplore.ieee.org/abstract/document/7285735
% 
% NOTE: the controller / WEC-Sim output is by default +F and +V into the actuator:
%                   _____________
% +F, +V --> ____________|       |
%                   _____|_______|
%
% but the HHEA (atleast this version) requires +F and +V out of the actuator
%  hence the negative sign when inputting controller.force and
%  controller.velocity into the hhea script
%                   _____________
% +F, +V <-- ____________|       |
%                   _____|_______|

clear;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load WEC-Sim cycle data using PI controller
% Output = structure containing:
% time (s)
% position, velocity, force - values of/on the actuator, in linear direction [m,m/s,N]
% ideal_velocity - experimental value related to the control theory. In rotational pitch direction [rad/s]
% applied_torque - torque applied by the actuator on the OSWEC flap [Nm]
% actual_angvel - The incoming angular velocity in pitch that is fed into the PI controller [rad/s]

% A = regular wave case
load ..\OSWEC\output_reg_nosum\wec_pi_out_A.mat
pi = output; clear output
load ..\OSWEC_Hydraulic_PTO\output_reg\wec_ptosim_out.mat
ptosim = output; clear output

% % B = irregular wave case
% load ..\OSWEC\output_irreg_nosum\wec_pi_out_B.mat
% pi = output; clear output
% load ..\OSWEC_Hydraulic_PTO\output_irreg\wec_ptosim_out.mat
% ptosim = output; clear output

% set time steps
pi.dt = pi.controller.time(2) - pi.controller.time(1);
ptosim.dt = ptosim.controller.time(2) - ptosim.controller.time(1);

% find power / work of the actuator.
% 'Ideal' in the sense that this is the maximum power input, and therefore 
% the max that the PTO can absorb
pi.idealPower = pi.controller.force .* pi.controller.velocity;
pi.idealWork = cumsum(pi.idealPower)*pi.dt;

ptosim.idealPower = ptosim.controller.force .* ptosim.controller.velocity;
ptosim.idealWork = cumsum(ptosim.idealPower)*ptosim.dt;

% calculate some mean quantities for PTO-Sim
ptosim.hydraulicMotor.torque = ptosim.ptosim.rotaryGenerator.genPower ./ ptosim.ptosim.hydraulicMotor.angVel;
ptosim.meanAbsPower = mean(ptosim.ptosim.pistonNCF.absPower(100/ptosim.dt:end));
ptosim.meanElecPower = mean(ptosim.ptosim.rotaryGenerator.elecPower(100/ptosim.dt:end));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the HHEA script using the PI controller input
hhea = struct();
[hhea.totalCost, hhea.M, hhea.I, hhea.genP, hhea.choices, hhea.scale, hhea.params] = ...
    HHEA_pto(-pi.controller.force,-pi.controller.velocity,pi.controller.time);

hhea.time = pi.controller.time;
hhea.meanAbsPower = mean(pi.idealPower(100/pi.dt:end)); % average actuator power after ramping period
hhea.meanElecPower = mean(hhea.genP(100/pi.dt:end)); % average elec power after ramping period
hhea.choices_label = {'DV_Ind_HECM1', 'PRA1', 'PRB1'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-process and plot 
set(groot, 'defaultLinelineWidth', 2.0);
set(0,'DefaultFigureWindowStyle','docked');

figure()
plot(pi.controller.time, pi.controller.force*1e-6,...
    pi.controller.time, pi.controller.velocity);
xlabel('Time (s)');
ylabel('Force (MN), velocity (rad/s)');
legend('Controller force (MN)', 'Velocity (rad/s)');
title('Actuator force and velocity required by PI controller');

figure()
p = plot(ptosim.controller.time, -ptosim.ptosim.pistonNCF.absPower*1e-3,...
     ptosim.controller.time, -ptosim.ptosim.rotaryGenerator.elecPower*1e-3, '--',...
     hhea.time, pi.idealPower*1e-3,...
     hhea.time, hhea.genP*1e-3, '--',...
     pi.controller.time, pi.wave.elevation,'k');
p(1).Color = [234/255, 112/255, 13/255];
p(2).Color = [234/255, 112/255, 13/255];
p(3).Color = [56/255, 112/255, 159/255];
p(4).Color = [56/255, 112/255, 159/255];
p(5).Color = [0/255, 0/255, 0/255];
legend('PTO-Sim Absorbed Power', 'PTO-Sim Electrical Power',...
    'HHEA Absorbed Power (PI)','HHEA Electrical Power','Wave elevation');
xlabel('Time (s)');
ylabel('Power (kW), Elevation (m)');
title('PTO-Sim and HHEA electrical power comparisons');

figure()
plot(pi.controller.time, pi.idealWork*1e-3,...
     hhea.time, cumsum(hhea.genP)*pi.dt*1e-3,...
     ptosim.controller.time, ptosim.idealWork*1e-3,...
     ptosim.controller.time, -cumsum(ptosim.ptosim.rotaryGenerator.elecPower)*ptosim.dt*1e-3);
title('Integral of work over time [kW]');
legend('HHEA actuator power (PI)','HHEA electrical power',...
    'PTO-Sim acuator power','PTO-Sim electrical output');

% print mean electrical output and efficiency over the entire cycle
fprintf('\nMean Electrical Power (kW): \nRM Power: %1.2f (%4.1f%s) \nHHEA Power: %1.2f, (%4.1f%s) \n',...
    ptosim.meanElecPower*1e-3,ptosim.meanElecPower/ptosim.meanAbsPower*100,'%',...
    hhea.meanElecPower*1e-3,hhea.meanElecPower/hhea.meanAbsPower*100,'%');

