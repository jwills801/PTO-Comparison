% This code is meant to maximize the power output from the check valve PTO.
% The assumption here is the the accumulators are large enough to absorbe
% any volume changes. The final volume is made equal to the ending volume
% only by use of the motor. The motor is sized to make these volumes be
% equal.

% The simulink model used for this code is OSWEC_Constant_DeltaP.slx
% This simulink model does not actually model the accumulator dynamics, it
% just provides some constant delta P, as dictated by this code

deltaP_vals = linspace(0,2e8,500); % Pa
Work_mech = NaN(size(deltaP_vals)); % initilize vector

for ind = 1:length(deltaP_vals)
    deltaP = deltaP_vals(ind);
    wecSim
    controller = struct();
    controller.velocity = controller1_out.signals.values(:,2); % + into actuator - piston velocity
    controller.force = controller1_out.signals.values(:,3); % + into actuator - piston force
    inst_power_mech = -controller.velocity.*controller.force;
    Work_mech(ind) = sum(inst_power_mech)*simu.time(2);
end

figure, plot(deltaP_vals,Work_mech/1e3), xlabel('Delta P'), ylabel('Mechanical Work In (kJ')

return
%% append controller output to controller input
controller = struct();
controller.time = simu.time;
controller.position = controller1_out.signals.values(:,1); % + towards actuator - piston position
controller.velocity = controller1_out.signals.values(:,2); % + into actuator - piston velocity
controller.force = controller1_out.signals.values(:,3); % + into actuator - piston force
controller.flow = controller1_out.signals.values(:,4); % flow into high Pressure Accum
controller.applied_torque = controller1_out.signals.values(:,5); % + towards actuator - torque on the WEC
controller.actual_angvel = controller1_out.signals.values(:,6); % + towards actuator - angular velocity of WEC
controller.signal_names = ["position" "velocity" "force" ... 
    "Flow" "WEC torque" "WEC angular velocity"];

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

inst_power_mech = -controller.velocity.*controller.force;
Work_mech = cumsum(inst_power_mech)*controller.time(2);
figure
plot(controller.time,Work_mech/1e3)
xlabel('Time (s)');
ylabel('Work (kJ)');
return

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






function Important_Quantities = Analyze(output)
    
    
    
    w = 3000*2*pi/60;
    inst_power_out = -.9*w*controller1_out.signals.values(:,6).*(controller1_out.signals.values(:,6)>=0) - 1/.9*w*controller1_out.signals.values(:,6).*(controller1_out.signals.values(:,6)<0) ;
    figure, plot(controller1_out.time,inst_power_out)

    figure, plot(controller1_out.time,cumsum(inst_power_out)*simu.dt)
    
    
    

end