% This code uses a genetic algorithm (built into matlab) to set parameters
% involved with the check valve PTO
% Specifically, these parameters are: 
%   1. volume of low pressure accum [m^3]
%   2. volume of high pressure accum [m^3]
%   3. pre-charge pressure of low pressure accum [Pa]
%   4. pre-charge pressure of high pressure accum [Pa]
%   5. displacement of motor [cc]

% The motor always spins at 3000 RPM. The electric efficiency is always 90%
% but the hydraulic efficiency varies with a pump map from pmouvahead.


nvars = 5; options = optimoptions('ga','PlotFcn', @gaplotbestf);
lb = [0;0;0;0;0]; % Lower bound on each design parameter
ub = [50; 50; 40e6; 40e6; 5000]; % Upper bound on each design parameter
[X,fval] = ga(@fun,nvars,[],[],[], [],lb,ub,[],options);


wecSim
%% append controller output to controller input
controller = struct();
controller.time = simu.time;
controller.position = controller1_out.signals.values(:,1); % + towards actuator
controller.velocity = controller1_out.signals.values(:,2); % + into actuator
controller.force = controller1_out.signals.values(:,3); % + into actuator
controller.applied_torque = controller1_out.signals.values(:,5); % + towards actuator
controller.actual_angvel = controller1_out.signals.values(:,6); % + towards actuator
controller.Q = controller1_out.signals.values(:,7);
controller.torque_generator = controller1_out.signals.values(:,8);
controller.signal_names = ["position" "velocity" "force" ... 
    "applied torque" "actual angular velocity"];
check = struct();
check.f = output.ptosim.pistonNCF.force;
% check.f2 = pto1_out.signals.values(:,);
check.signal_names = ["position" "velocity" "acceleration" "forceTotal" ...
    "forceActuation" "forceConstraint" "forceInternalMechanics" "powerInternalMechanics"];
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

inst_power = -controller.velocity.*controller.force;
inst_power_out = -.9*w*controller.torque_generator.*(controller.torque_generator>=0) - 1/.9*w*controller.torque_generator.*(controller.torque_generator<0) ;
t = controller.time;
Work_in = cumsum(inst_power)*controller.time(2);
Work_out = cumsum(inst_power_out)*controller.time(2);
figure
plot(t,Work_in/1e6,t,Work_out/1e6)
xlabel('Time (s)');
ylabel('Work (MJ)');
legend('Work In','Work Out');

figure
plot(controller.time,inst_power_out/1000,controller.time,inst_power/1000)
legend('Power out','Power in')
xlabel('Time (s)');
ylabel('Power (kW)');
xlim([100 125])

%% Find size of Hydraulic pump motor
MotorSize = d_scale*107 %cc
genSize = max(abs(inst_power_out))/1e3 % kW
Efficiency = sum(inst_power_out)/sum(inst_power)

B = 18; % m (This is the width of the oswec)
t_start = 50;
dt= output.ptosim.time(2);
t= output.ptosim.time;
Energy_in_first_chunk = -sum(inst_power(1:find(t==t_start)))*dt;
Energy_out_first_chunk = -sum(inst_power_out(1:find(t==t_start)))*dt;
Total_energy_in = sum(inst_power)*dt;
Total_energy_out = sum(inst_power_out)*dt;
ave_power_in = (Total_energy_in - Energy_in_first_chunk)/150;
ave_power_out = (Total_energy_out - Energy_out_first_chunk)/150;
CWR_in = ave_power_in/waves.Pw/B
CWR_out = ave_power_out/waves.Pw/B


%% plot accumulator dynamics
figure, plot(output.ptosim.time,output.ptosim.accumulator(1).pressure,output.ptosim.time,output.ptosim.accumulator(2).pressure)
legend('High Pressure','Low Pressure')
xlabel('Time (s)')
ylabel('Pressure (Pa)')

figure, plot(output.ptosim.time,output.ptosim.accumulator(1).volume,output.ptosim.time,output.ptosim.accumulator(2).volume)
legend('High Pressure','Low Pressure')
xlabel('Time (s)')
ylabel('Change in Volume (m^3)')





return

%%
X = [5 5 1e6 10e6 1200];
fun(X)

function J = fun(X)
    wecSim
    Volume_punishment = abs(output.ptosim.accumulator(1).volume(end)) ; % punish not ending the simulation with the same volume as it began with
    
    inst_power_out = -.9*w*controller1_out.signals.values(:,8).*(controller1_out.signals.values(:,8)>=0) - 1/.9*w*controller1_out.signals.values(:,8).*(controller1_out.signals.values(:,8)<0) ;
    Energy_reward = sum(inst_power_out)*simu.dt;
    Energy_reward = Energy_reward/.5e8 ; % scaled so that maximum overall reward is around 2
    J = Energy_reward - Volume_punishment;
    

end