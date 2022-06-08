% Script used to compare the power output of the RM5 PTO and the HHEA
% system. 
% Utilizes the output of a controller used as a PTO in wec-sim. Controller
% must output force on and velocity of an actuator. 
% Take controller-wecSim output, RMP-wecSim output and hhea(f,v)

clc;
close all;
clear;

tic

% Load WEC-Sim cycle data using PI controller
% _A for regular case, _B for irregular
load ../wec-PI_out_A_alternate.mat;
dt = output.controller.time(2) - output.controller.time(1);
simu.dtOut = dt;
simu.dt = dt;
idealPower = output.controller.force .* -output.controller.velocity;
idealE = cumsum(idealPower)*dt;

% Load RMP data from a full wec-sim run.
load rmp_equivalent_A.mat
rmp.time = output.controller.time;
rmp.hydMot.torque = rmp.rotGen.genPower ./ rmp.hydMot.angVel;
rmp.meanAbsPower = mean(rmp.piston.absPower(100/dt:end));
rmp.meanElecPower = mean(rmp.rotGen.elecPower(100/dt:end));

toc

hhea = struct();
hhea.time = output.controller.time;

%Jackson Adjusted this line, as well as the function HHEA_pto
[hhea.totalCost, hhea.M, hhea.I, hhea.genP, hhea.choices, hhea.scale, hhea.params] = ...
    HHEA_pto(output.controller.force,-output.controller.velocity,output.controller.time);
hhea.meanAbsPower = mean(idealPower(100/dt:end)); % average actuator power after ramping period
hhea.meanElecPower = mean(hhea.genP(100/dt:end)); % average elec power after ramping period
hhea.choices_label = {'DV_Ind_HECM1', 'PRA1', 'PRB1'};

choices = struct;
choices.time = hhea.time;
choices.signals.values = hhea.choices;
choices.signals.dimensions = [length(hhea.choices), 3];
choices.label = hhea.choices_label;
choices.dt = dt;
save('hhea_choices.mat', 'choices');

toc





% notes: 
% - calculate efficiency based on MP loss after optimization -> done
% - don't calculate total power output based on losses (may be
% approximations) -> done
% - check if DCV is even being used

%% Jackson's Additions
t = output.controller.time;
figure(10)
plot(t, hhea.battery)

%% Analysis
set(groot, 'defaultLinelineWidth', 2.0);
set(0,'DefaultFigureWindowStyle','docked');

figure()
plot(output.controller.time, output.controller.force*1e-6,...
    output.controller.time, output.controller.velocity);
xlabel('Time (s)');
ylabel('Force (MN), velocity (rad/s)');
legend('Controller force (MN)', 'Velocity (rad/s)');
title('Actuator force and velocity required by PI controller');

figure()
plot(rmp.time, rmp.hydMot.angVel, rmp.time, rmp.hydMot.torque);
legend('ang speed', 'torque');
title('RMP motor');


% cycle = load('../JCB_LongerDC.mat');
% figure()
% plot(hhea.time, -output.controller.force,...
%     cycle.Grading.time, cycle.Grading.boom_f);
% legend('WEC force','JCB force');
% xlabel('Time (s)');
% ylabel('Force (N)');
% title('wec and jcb cycle comparison - force');

% figure()
% plot(hhea.time, -output.controller.velocity*hhea.scale,...
%     cycle.Grading.time, cycle.Grading.boom_v);
% legend('WEC velocity*hhea.scale', 'JCB velocity');
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% title('wec and jcb cycle comparison - velocity');


% figure()
% plot(hhea.time, idealPower, hhea.time, hhea.M, '--');
% xlabel('Time (s)');
% ylabel('Power (W)');
% legend('Actuator input power', 'HHEA Losses');
% title("HHEA input power and losses. Scale = " + num2str(hhea.scale));

% figure()
% plot(rmp.time, rmp.piston.absPower, 'b',...
%      rmp.time, rmp.rotGen.elecPower, 'b--',...
%      rmp.time(100/dt:end), rmp.meanElecPower*ones(size(rmp.time(100/dt:end))),'b',...
%      rmp.time, idealPower, 'c',...
%      hhea.time, hhea.power, 'c--',...
%      hhea.time(100/dt:end), hhea.meanElecPower*ones(size(hhea.time(100/dt:end))),'c');
% legend('RM5 Abs. Power', 'RM5 PTO Elec Power','RM5 post-ramp mean power',...
%     'HHEA Ideal Power (PI)','HHEA Elec Power','HHEA post-ramp mean power');
% xlabel('Time (s)');
% xlim([0 55]);
% ylabel('Power (W)');
% title('Ideal, RM5, and HHEA electrical power comparisons');


% figure()
% p = plot(rmp.time, rmp.piston.absPower*1e-3,...
%      rmp.time, rmp.rotGen.elecPower*1e-3, '--',...
%      hhea.time, idealPower*1e-3,...
%      hhea.time, hhea.genP*1e-3, '--',...
%      hhea.time, wave_elevation);
% p(1).Color = [234/255, 112/255, 13/255];
% p(2).Color = [234/255, 112/255, 13/255];
% p(3).Color = [56/255, 112/255, 159/255];
% p(4).Color = [56/255, 112/255, 159/255];
% p(5).Color = [0/255, 0/255, 0/255];
% legend('RM5 Absorbed Power', 'RM5 Electrical Power',...
%     'HHEA Absorbed Power (PI)','HHEA Electrical Power','Wave elevation * 100');
% xlabel('Time (s)');
% xlim([30 81.5]);
% ylabel('Power (kW), Elevation (m)');
% % title('Ideal, RM5, and HHEA electrical power comparisons');
% 
% set(gca, 'FontSize', 14);
% set(gca, 'FontWeight', 'Bold');

figure()
plot(hhea.time,hhea.params.qr1,...
    hhea.time,hhea.params.qr2,...
    hhea.time,hhea.params.qr3,...
    hhea.time,hhea.params.hecmFlow,...
    hhea.time,hhea.params.mpFlow);
legend('pressure rail 1','pressure rail 2','pressure rail 3',...
    'hecm','main pump');
title('flow rates');

figure()
plot(hhea.time,hhea.params.qr1,...
    hhea.time,hhea.params.qr2,...
    hhea.time,hhea.params.qr3);
legend('pressure rail 1','pressure rail 2','pressure rail 3');
title('flow rates');

figure()
plot(hhea.time, cumsum(idealPower)*1e-3,...
     hhea.time, cumsum(hhea.genP)*1e-3);
title('integral of work over time');
legend('input power','electrical output');

fprintf('\nMean Electrical Power (kW): \nRM Power: %1.2f (%4.1f%s) \nHHEA Power: %1.2f, (%4.1f%s) \n',...
    rmp.meanElecPower*1e-3,rmp.meanElecPower/rmp.meanAbsPower*100,'%',...
    hhea.meanElecPower*1e-3,hhea.meanElecPower/hhea.meanAbsPower*100,'%');
fprintf('\nStd HHEA cycle: Qmax = %1.5f, Pmax = %4.1f \n',...
    max(abs(cycle.Grading.boom_v))*cycle.boom_Ac,...
    max(abs(cycle.Grading.boom_f))/cycle.boom_Ar*1e-6);
fprintf('RMP/wec-sim cycle:  Qmax = %1.5f, Pmax = %4.1f \n\n\n',...
    max(abs(output.controller.velocity))*capA,...
    max(abs(rmp.acc(1).pressure-rmp.acc(2).pressure))/capA*1e-6);

