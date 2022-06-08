%% EHA Case
%% Define Force and velocity from WEC-sim
load('wec-PI_out_A_alternate.mat')
F = output.controller.force;
v = output.controller.velocity;
t = output.controller.time;
dt = diff(t);
A = 0.0378;

%% Calculate Actual Flow and delta P through/across variable displacement motor

%Calculate Flow
Q_Act = v*A;
n = length(Q_Act);

% Pressure Differential (Pressure actuator - Pressure Rail)
deltaP = F/A;


%% Fluid Properties

mu=(32e-6)*870;
B = 1.7e9;
rho = 870;

%% Define Pump Constants

%% I am scaling the pump to be X times larger than the 107 cc b/c That will make -1<fracDisp<1
Scale = 6;  % Regular wave case
%Scale = 20; % for irregular wave case
%% Manufacturer 107cc/rev
% Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
    D = 107; % cc/rev
    d = (D*100^-3)/(2*pi); % m^3/rad 
% Torques Loss Constants
    Cf =  53.7e-3;
    Ch = 53.6;
    Cv = 23.5e3;
% Flow Loss Constants
    Cs = 4.26e-9;
    Cst = 0*1e-5;

%% Define Pump Parameters to Study
% (fractional displacement/angular velocity/pressure differential)

% Angular Velocity
Wrpm = 3000; %revolutions per minute
w = Wrpm.*(2*pi/60); % radians per second



%% Find fraction of total displacement
fracDisp = NaN(n,1);
T_Ideal = NaN(n,1);
TLoss = NaN(n,1);
T_Act = NaN(n,1);
P_out = NaN(n,1);
P_in = NaN(n,1);
QLoss = NaN(n,1);
P_leak = NaN(n,1);
P_mech = NaN(n,1);
P_L_elect = NaN(n,1);
P_out_sum = 0;
P_in_sum = 0;



% Need to check what quadrant this is in, because the equations differ
% QLoss and TLoss scale with the pump. Pump 10X larger => Qloss 10X larger 
for i = 1:n
    % Calculate fracDisp assuming it is a posative value
    fracDisp(i) = (Q_Act(i) + sign(deltaP(i))*Scale*abs(d*Cs*(deltaP(i))/mu) + sign(deltaP(i))*Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5))/(w*d*Scale-sign(deltaP(i))*Scale*abs(d*w*deltaP(i)/B));
        if fracDisp(i) <= 0 % if the assumption that fracdisp is + is incorrect, recalculate fracDisp assuming fracDisp is - 
            fracDisp(i) = (Q_Act(i) + sign(deltaP(i))*Scale*abs(d*Cs*(deltaP(i))/mu) + sign(deltaP(i))*Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5))/(w*d*Scale+sign(deltaP(i))*Scale*abs(d*w*deltaP(i)/B));
        end
    
    % Check if fracDisp was calculated correctly
    % Calculate the error due to the fractional displacement calculation
    % If I use the fracDisp I calculated to solve for Q_Act (which was given)
    % Do I get the same value?
    QLoss(i) = Scale*abs(d*Cs*(deltaP(i))/mu) + Scale*abs(fracDisp(i)*d*w*(deltaP(i))/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5);
    Q_Act_calc(i) = w*d*fracDisp(i)*Scale - sign(deltaP(i))*QLoss(i);
    
    
    T_Ideal(i) = deltaP(i)*d*fracDisp(i)*Scale;
    TLoss(i) = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP(i))*Cf) + abs(fracDisp(i)*Ch*w^2*rho*d^(5/3)/2)  );
    
    
    
    Q_Ideal(i) = fracDisp(i)*d*w*Scale;
    T_Act(i) = T_Ideal(i) + sign(w)*TLoss(i); % this is a plus sign because if T_Ideal is always neg, then the TLoss needs to decrease its magnitude
    
    % Power out with 90% effiency
    if T_Act(i) < 0
        P_out(i) = .9*w*T_Act(i);
        P_L_elect(i) = -.1*w*T_Act(i); % the negative sign makes sure the loss is posative 
    else
        P_out(i) = w*T_Act(i)/.9;
        P_L_elect(i) =  .111*w*T_Act(i); % This loss accounts for energy that needs to come FROM the generator to the system.
    end
    
    P_out_sum = P_out_sum + P_out(i);
    P_in(i) = (F(i)*v(i)); 
    P_in_sum = P_in(i) + P_in_sum;
    
    
    % check if the difference in work in and work out  (i.e. integral of the power between in and out) is the losses in the hydraulic and electrical components.
    P_leak(i) = abs(deltaP(i) *QLoss(i));
    P_mech(i) = w*TLoss(i);
   
end

e = P_out_sum/P_in_sum;
EHA_Generator = max(abs(P_out)); % W
EHA_Pump_size = D*Scale;

figure
hold on
scatter(fracDisp,deltaP,'r','x')
run('../NEWresults_one/EfficiencyMap.m')
hold off

%return
figure(2)
plot(t,-P_out,t,-P_in)
legend('Power Out','Power In')
xlabel('Time (s)')
ylabel('Power (J)')

dt = mean(diff(t));
figure(1)
plot(t,-[cumsum(P_in),cumsum(P_out)]);
legend('Work in', 'Work out');
xlabel('Time(s)');
ylabel('Work (J)');
title(['Efficiency = ',num2str(sum(P_out)/sum(P_in)*100),'%'])
mean_power_out = -P_out_sum/length(t);
max_power_out = max(-P_out);
mean_power_in = -P_in_sum/length(t);
max_power_in = max(-P_in);
%% Regular Wave Case
% mean power out is 327kW; max power out is 836kW (256% sized) 
% mean power in is 430kW; max power in is 1060kW (246% sized)   
% Pump/motor is about 107*11=1177cc 
% efficiency = 74.25%

%% Irregular Wave Case
% mean power out is 196kW; max power out is 2399kW (1227% sized) 
% mean power in is 301kW; max power in is 3091kW (1028% sized)   
% Pump/motor is about 107*20=2140cc 
% efficiency = 65.04%

% check if fracDisp was calculated correctly
figure(3)
plot(t,Q_Act,t,Q_Act_calc)
legend('Q_Act','calculated with fracDisp')
xlabel('Time (s)')
ylabel('Flow (m^3/s)')

% Check if the individual sources of power loss add up to the total power loss?
figure(4)
plot(t,(P_in-P_out),t,-(P_leak+P_mech+P_L_elect))
legend('Power loss','Loss due to Flow,torque and effiency')
xlabel('Time (s)')
ylabel('Power (J)')


figure(5)
plot(t,P_leak,t,P_mech,t,P_L_elect,t,-(P_in-P_out))
legend('Leakage','friction','efficiency of generator')
xlabel('Time (s)')
ylabel('Power (J)')
title('Power Loss Contributions over time')


figure(7)
plot(t,P_leak,t,P_mech)
legend('leak','mech')
title('Power Loss Contributions over time')

figure(8)
plot(t,fracDisp)
% title('1177cc Pump Motor')
xlabel('Time (s)')
ylabel('Fractional Displacement')
grid on
print('fracdisp', '-dpng', '-r600')

figure(9)
subplot(2,1,1);
plot(t,-F)
title('Force Trajectory')
ylabel('Force (N)')
xlabel('Time (s)')
grid on
subplot(2,1,2); 
plot(t,v)
title('Velocity Trajectory')
ylabel('Velocity (m/s)')
xlabel('Time (s)')
grid on
print('DriveCycle', '-dpng', '-r600')

figure(10)
plot(t,-cumsum(F.*v)*dt(1))
title('Energy into system')
ylabel('Work (J)')
xlabel('Time (s)')

figure(11)
yyaxis left
plot(t,-P_in,t,-P_out)
xlabel('Time (s)')
ylabel('Power (W)')
ylim([-2e5 12e5])
yyaxis right
plot(t,v)
legend('Power In','Power Out','Velocity','location','northwest')
xlim([150 200])
ylim([-2 12])
ylabel('Velocity (m/s)')
grid on
print('EHAPower', '-dpng', '-r600')