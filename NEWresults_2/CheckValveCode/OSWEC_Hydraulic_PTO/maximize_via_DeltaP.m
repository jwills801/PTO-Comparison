% This code is meant to maximize the power output from the check valve PTO.
% The assumption here is the the accumulators are large enough to absorbe
% any volume changes. The final volume is made equal to the ending volume
% only by use of the motor. The motor is sized to make these volumes be
% equal.

% The simulink model used for this code is OSWEC_Constant_DeltaP.slx
% This simulink model does not actually model the accumulator dynamics, it
% just provides some constant delta P, as dictated by this code

% For regular Waves, most energy is aborbed using deltaP of 8.4e6 Pa.
% Mechanical Energy absorbed: 98 MJ
% GenSize: 509 kW
% Motor Size: 1554 cc
% CWR_in: .7081
% CWR_out: .5070

% For Irregular Waves, most energy is aborbed using deltaP of 8.1e6 Pa.
% Mechanical Energy absorbed: 44.4 MJ
% GenSize: 231 kW
% Motor Size: 731 cc
% CWR_in: 0.8843
% CWR_out: 0.5485

deltaP_vals = linspace(0,2e7,100); % Pa
deltaP_vals = 8.1e6; % Pa
Work_mech = NaN(size(deltaP_vals)); % initilize vector

outer_toc = tic;
for ind = 1:length(deltaP_vals)
    deltaP = deltaP_vals(ind);
    evalc('wecSim');
    controller = struct();
    controller.velocity = controller1_out.signals.values(:,2); % + into actuator - piston velocity
    controller.force = controller1_out.signals.values(:,3); % + into actuator - piston force
    inst_power_mech = -controller.velocity.*controller.force;
    Work_mech(ind) = sum(inst_power_mech)*simu.time(2);
    
    disp([num2str(ind) ' of ' num2str(length(deltaP_vals)) ' simulations complete'])
end
toc(outer_toc)

[Best_Energy_Capture, ind]= max(Work_mech);
Best_DeltaP = deltaP_vals(ind);
close all
figure, plot(deltaP_vals,Work_mech/1e6), xlabel('Delta P (Pa)'), ylabel('Mechanical Work In (MJ)')

%% Run wecsim again with optimal deltaP in order to get the whole picture
deltaP = Best_DeltaP;
evalc('wecSim');
controller = struct();
controller.time = simu.time;
controller.position = controller1_out.signals.values(:,1); % + towards actuator - piston position
controller.velocity = controller1_out.signals.values(:,2); % + into actuator - piston velocity
controller.force = controller1_out.signals.values(:,3); % + into actuator - piston force
controller.flow = abs(controller1_out.signals.values(:,4)); % flow into high Pressure Accum
controller.applied_torque = controller1_out.signals.values(:,5); % + towards actuator - torque on the WEC
controller.actual_angvel = controller1_out.signals.values(:,6); % + towards actuator - angular velocity of WEC

T = controller.time(end);
T = T - simu.rampTime; % How much time does the main pump have to pump?
Net_Flow = sum(controller.flow)*controller.time(2); % integrate flow to get volume (m^3)
Ave_Flow = Net_Flow/T; % Motor can remove net flow (m^3/s)
w = 2000*2*pi/60; % Angular Velocity of e generator [rad/s]
D = Ave_Flow/w*2*pi*1e6;  % Size of motor (cc)
disp(['Pump size is ', num2str(D,4),' cc'])

Constant_eff = .78; % From Pmouvahead at 35 MPa and 2000 RMP and full displacement

%% Plots
%Plot waves
waves.plotEta(simu.rampTime);
try 
    waves.plotSpectrum();
catch
end
figure
subplot(2,1,1);
plot(controller.time,controller.force)
ylabel('Force (N)')
xlabel('Time (s)')
grid on
subplot(2,1,2); 
plot(controller.time,controller.velocity)
ylabel('Velocity (m/s)')
xlabel('Time (s)')
grid on

t= controller.time;
dt = controller.time(2);
inst_power_mech = -controller.velocity.*controller.force;
Work_mech = cumsum(inst_power_mech)*dt;
inst_power_elec = [zeros(simu.rampTime/dt+1,1);ones(T/dt,1)]*Constant_eff*Ave_Flow * deltaP; %eta * flow * pressure
Work_elec = cumsum(inst_power_elec)*dt;
figure
plot(t,Work_mech/1e6,t,Work_elec/1e6)
xlabel('Time (s)');
ylabel('Work (MJ)');
legend('Work in','Work out')

figure
plot(controller.time,inst_power_mech/1000,controller.time,inst_power_elec/1000)
legend('Power in','Power out')
xlabel('Time (s)');
ylabel('Power (kW)');
xlim([100 125])

%% Find size of Hydraulic pump motor
genSize = max(abs(inst_power_elec))/1e3 % kW

B = 18; % m (This is the width of the oswec)
t_start = 50;
t= controller.time;
dt = t(2);
Energy_in_first_chunk = sum(inst_power_mech(1:find(t==t_start)))*dt;
Energy_out_first_chunk = sum(inst_power_elec(1:find(t==t_start)))*dt;
Total_energy_in = sum(inst_power_mech)*dt;
Total_energy_out = sum(inst_power_elec)*dt;
ave_power_in = (Total_energy_in - Energy_in_first_chunk)/(t(end)-t_start);
ave_power_out = (Total_energy_out - Energy_out_first_chunk)/(t(end)-t_start);
CWR_in = ave_power_in/waves.Pw/B
CWR_out = ave_power_out/waves.Pw/B


%% Plot excitation force
a =5; b = 3.9;
x = NaN(length(t),2);
x(:,1) = output.bodies(1).position(:,5); % Initial position in pitch
x(:,2) = output.bodies(1).velocity(:,5); % Initial velocity in pitch
moment_arm = [a*cos(x(:,1)) zeros(size(x(:,1))) -a*sin(x(:,1)) zeros(size(x(:,1))) ones(size(x(:,1))) zeros(size(x(:,1)))]; % Converts forces in 6 DOF to the only DOF we care about: Pitch
excitation = sum(moment_arm.*output.bodies(1).forceExcitation(:,:),2);

figure, plot(t,excitation/1e6), ylabel('Excitation Torque [MNm]'), xlabel('Time [s]'), xlim([100 300]), grid on
