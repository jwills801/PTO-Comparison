%% Housekeeping
return
%close all
%clear table  

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

% Plot RY forces for body 1
%plotForces(output,1,5)

%Plot RY response for body 1
%output.plotResponse(1,5);

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
%figure, plot(output.ptosim.hydraulicMotor(1).angVel,output.ptosim.hydraulicMotor(1).volFlowM), xlabel('angular velocity'), ylabel('Flow')
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



%% Plots 
% set(0,'DefaultFigureWindowStyle','docked')
% 
% figure();
% plot(output.ptosim.time,output.ptosim.pistonNCF.topPressure/1e6,output.ptosim.time,output.ptosim.pistonNCF.bottomPressure/1e6)
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Pressure (MPa)')
% title('Piston Pressure')
% legend('topPressure','bottomPressure')
% grid on
% 
% figure();
% plot(output.ptosim.time,output.ptosim.pistonNCF.force/1e6)
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Force (MN)')
% title('PTO Force')
% grid on
% 
% 
% figure();
% plot(output.ptosim.time,output.ptosim.pistonNCF.absPower/1e3,output.ptosim.time,output.ptosim.rotaryGenerator.genPower/1e3,output.ptosim.time,output.ptosim.rotaryGenerator.elecPower/1e3)
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Power (kW)')
% title('Absorbed Power, Mechanical Power, and Electrical Power')
% legend('absPower','mechPower','elecPower')
% grid on
% 
% figure();
% plot(output.ptosim.time,output.ptosim.hydraulicMotor.angVel)
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Speed (rad/s)')
% title('Speed')
% grid on

%% figure()
% subplot(3,1,1)
% plot(output.ptosim.time,(output.ptosim.accumulator(1).pressure-output.ptosim.accumulator(2).pressure)/1e6)
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Pressure (MPa)')
% title('Pressure Differential Between the Two Accumulators')
% grid on
% subplot(3,1,2)
% plot(output.ptosim.time,output.ptosim.hydraulicMotor.volFlowM./output.ptosim.hydraulicMotor.angVel)
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Volume (m^3)')
% title('Hydraulic Motor Volume')
% grid on
% subplot(3,1,3)
% plot(output.ptosim.time,((output.ptosim.accumulator(1).pressure-output.ptosim.accumulator(2).pressure).*(output.ptosim.hydraulicMotor.volFlowM./output.ptosim.hydraulicMotor.angVel))/(ptosim.rotaryGenerator.desiredSpeed*1.05))
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Damping (kg-m^2/s)')
% title('Generator Damping')
% grid on