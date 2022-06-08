clear
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
Scale = 13;  % Regular wave case
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
Wrpm = 1800; %revolutions per minute
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
    %if deltaP(i) >= 0
        fracDisp(i) =   ( Q_Act(i) - Scale*abs(d*Cs*(deltaP(i))/mu) - Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5)  )/(w*d*Scale+Scale*abs(d*w*deltaP(i)/B));
        if fracDisp(i) <= 0
            fracDisp(i) = (Q_Act(i) - Scale*abs(d*Cs*(deltaP(i))/mu) - Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5))/(w*d*Scale-Scale*abs(d*w*deltaP(i)/B));
        end
    %else
       % fracDisp(i) = (Q_Act(i) - Scale*abs(d*Cs*(deltaP(i))/mu) - Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5))/(w*d*Scale+Scale*abs(d*w*deltaP(i)/B));
      %  if fracDisp(i) <= 0
      %      fracDisp(i) = (Q_Act(i) - Scale*abs(d*Cs*(deltaP(i))/mu) - Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5))/(w*d*Scale-Scale*abs(d*w*deltaP(i)/B));
      %  end
    %end
    
    QLoss(i) = Scale*abs(d*Cs*(deltaP(i))/mu) + Scale*abs(fracDisp(i)*d*w*(deltaP(i))/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5);

    
    T_Ideal(i) = deltaP(i)*d*fracDisp(i)*Scale;
    TLoss(i) = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP(i))*Cf) + abs(fracDisp(i)*Ch*w^2*rho*d^(5/3)/2)  );
    
    
    
    Q_Ideal(i) = fracDisp(i)*d*w*Scale;
    T_Act(i) = T_Ideal(i) + sign(w)*TLoss(i);
    
    % Power out with 90% effiency
    if T_Act(i) < 0
        P_out(i) = .9*w*T_Act(i);
        P_L_elect(i) = -.1*w*T_Act(i); % the negative sign makes sure the loss is posative 
    else
        P_out(i) = w*T_Act(i)/.9;
        P_L_elect(i) =  .111*w*T_Act(i);
    end
    
    P_out_sum = P_out_sum + P_out(i);
    P_in(i) = (F(i)*v(i)); 
    P_in_sum = P_in(i) + P_in_sum;
    
    
    % check if the difference in work in and work out  (i.e. integral of the power between in and out) is the losses in the hydraulic and electrical components.
    P_leak(i) = abs(deltaP(i) *QLoss(i));
    P_mech(i) = w*TLoss(i);
   
end

e = P_out_sum/P_in_sum
figure(2)
plot(t,-P_out,t,-P_in)
legend('Power Out','Power In')
xlabel('Time (s)')
ylabel('Power (J)')

dt = mean(diff(t));
figure(1)
plot(t,[cumsum(P_in),cumsum(P_out)]*dt);
legend('Work in', 'Work out');
xlabel('Time(s)');
ylabel('Cumulative energy in and out (J)');
title(['Efficiency = ',num2str(sum(P_out)/sum(P_in)*100),'%'])
mean_power_out = -P_out_sum/length(t);
max_power_out = max(-P_out);
mean_power_in = -P_in_sum/length(t);
max_power_in = max(-P_in);
%% Regular Wave Case
% mean power out is 338kW; max power out is 918kW (272% sized) 
% mean power in is 430kW; max power in is 1060kW (246% sized)   
% Pump/motor is about 107*13=1391cc 
% efficiency = 78.61%

%% Irregular Wave Case
% mean power out is 216kW; max power out is 2512kW (1161% sized) 
% mean power in is 301kW; max power in is 3091kW (1028% sized)   
% Pump/motor is about 107*20=2140cc 
% efficiency = 71.99%


%% CHECKS
%% Check if fracDisp was calculated correctly
% Calculate the error due to the fractional displacement calculation
% If I use the fracDisp I calculated to solve for Q_Act (which was given)
% Do I get the same value?
error = NaN(n,1);
Q_Act_calc = NaN(n,1);
for j=1:n
    if deltaP(i) >= 0
        Q_Act_calc(j) = w*d*fracDisp(j)*Scale - QLoss(j);
    else
        Q_Act_calc(j) = w*d*fracDisp(j)*Scale + QLoss(j);
    end
    error(j) = Q_Act(j)-Q_Act_calc(j);
end
figure(3)
plot(t,Q_Act,t,Q_Act_calc)
legend('Q_Act','calculated with fracDisp')
xlabel('Time (s)')
ylabel('Flow (m^3/s)')

%% Check if the individual sources of power loss add up to the total power loss?
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

figure(6)
plot(t,(T_Act*w-P_in),t,P_mech+P_leak); legend('Difference','sum')


figure(7)
plot(t,P_leak,t,P_mech)
legend('leak','mech')