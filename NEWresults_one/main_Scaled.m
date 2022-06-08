%% EHA Case

%% Define Force and velocity from WEC-sim
load Li_PI_one.mat
F = myoutput.signals.values(:,17);
v = myoutput.signals.values(:,11);
t = myoutput.time;
dt = t(2) - t(1);
A = 0.2382; A = .2;


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
Scale = 45;  % Regular wave case
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
    % Do I get the same values?
    QLoss(i) = Scale*abs(d*Cs*(deltaP(i))/mu) + Scale*abs(fracDisp(i)*d*w*(deltaP(i))/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP(i))/rho)^.5);
    Q_Act_calc(i) = w*d*fracDisp(i)*Scale - sign(deltaP(i))*QLoss(i);
    
    
    T_Ideal(i) = deltaP(i)*d*fracDisp(i)*Scale;
    TLoss(i) = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP(i))*Cf) + abs(fracDisp(i)*Ch*w^2*rho*d^(5/3)/2)  );
    
    
    
    Q_Ideal(i) = fracDisp(i)*d*w*Scale;
    T_Act(i) = T_Ideal(i) + sign(w)*TLoss(i); % |T_Act| needs be < |T_Ideal|
    
    % Power out with 90% effiency
    if T_Act(i) < 0
        P_out(i) = .9*w*T_Act(i);
        P_L_elect(i) = -.1*w*T_Act(i); % the negative sign makes sure the loss is posative 
    else
        P_out(i) = w*T_Act(i)/.9;
        P_L_elect(i) =  (1/.9-1)*w*T_Act(i); % This loss accounts for energy that needs to come FROM the generator to the system.
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

% figure
% hold on
% scatter(fracDisp,deltaP,'r','x')
% EfficiencyMap
% hold off
%  
%%
close all
figure()
plot(t,-P_out,t,-P_in)
hold on
%legend('Power Out','Power In')
xlabel('Time (s)')
ylabel('Power (W)')
grid on
xlim([150 164])
x = [.33 .39]; y = [.85 .88]; annotation('textarrow',x,y,'String','Power In');
x = [.43 .48]; y = [.2 .24]; annotation('textarrow',x,y,'String','Power Out');
line([150 164],[0 0 ],'LineWidth',3,'Color','k')
hold off



return
posPower = (-F.*v/1000>=0).*(-F.*v/1000);
negPower = (-F.*v/1000<=0).*(-F.*v/1000);
figure
area(t,posPower)
hold on
area(t,negPower)
xlim([150 164])
xlabel('Time (s)')
ylabel('Power (kW)')
set(gca,'box','off') 
grid on
hold off



dt = mean(diff(t));
figure()
plot(t,-[cumsum(P_in)*dt,cumsum(P_out)*dt]);
legend('Work in', 'Work out');
xlabel('Time(s)');
ylabel('Work (J)');
title(['Efficiency = ',num2str(sum(P_out)/sum(P_in)*100),'%'])
mean_power_out = -P_out_sum/length(t);
max_power_out = max(-P_out);
mean_power_in = -P_in_sum/length(t);
max_power_in = max(-P_in);


% check if fracDisp was calculated correctly
% figure(3)
% plot(t,Q_Act,t,Q_Act_calc)
% legend('Q_Act','calculated with fracDisp')
% xlabel('Time (s)')
% ylabel('Flow (m^3/s)')

% Check if the individual sources of power loss add up to the total power loss?
figure()
plot(t,-(P_in-P_out),t,(P_leak+P_mech+P_L_elect))
legend('Power loss','Loss due to Flow,torque and effiency')
xlabel('Time (s)')
ylabel('Power (J)')


figure()
plot(t,fracDisp)
% title('1177cc Pump Motor')
xlabel('Time (s)')
ylabel('Fractional Displacement')
grid on
print('Figures/fracdisp', '-dpng', '-r600')

figure
subplot(2,1,1);
plot(t,-F)
title('Force Trajectory')
ylabel('Force (N)')
xlabel('Time (s)')
grid on
xlim([150 164])
subplot(2,1,2); 
plot(t,v)
title('Velocity Trajectory')
ylabel('Velocity (m/s)')
xlabel('Time (s)')
ylim([-1.1 1.1])
grid on
xlim([150 164])

figure
yyaxis left
plot(t,-P_in/1000,t,-P_out/1000)
xlabel('Time (s)')
ylabel('Power (kW)')
%ylim([-2e5 12e5])
yyaxis right
plot(t,v)
legend('Power In','Power Out','Velocity','location','northwest')
xlim([150 200])
ylim([-2 12])
ylabel('Velocity (m/s)')
grid on

figure()
plot(t,T_Act,t,T_Ideal)
ylabel('Torque (Nm)'), xlabel('Time (s)')
legend('Actual Torque','Ideal Torque')

figure()
plot(t,Q_Act,t,Q_Ideal)
ylabel('(m^3/s)'), xlabel('Time (s)')
legend('Actual flow','Ideal flow')

%figure(), plot(t,T_Act,t,deltaP*d*Scale.*fracDisp + sign(w).*TLoss), legend('True','sum'), ylabel('Torque'), title('Torque')

%figure(), plot(t,Q_Act,t,w*d*Scale.*fracDisp - sign(deltaP).*QLoss), legend('true','sum'), ylabel('Flow'), title('Flow')

figure, plot(t,cumsum(-P_in)*dt)

%% Instantaneous efficiency
Hydraulic_in = (-P_in>=0).*(-P_in);
Hydraulic_out = -(-P_in<=0).*(-P_in);
Mechanical_in = -(-P_out<=0).*(-P_out);
Mechanical_out = (-P_out>=0).*(-P_out);
inst_e = -(Hydraulic_out-Mechanical_out)./(Hydraulic_in-Mechanical_in);
figure, plot(t,inst_e), ylabel('Instantaneous efficiency'), xlabel('Time (s)')

ave_Hydraulic_in = sum((-P_in>=0).*(-P_in))*dt;
ave_Hydraulic_out = sum(-(-P_in<=0).*(-P_in))*dt;
ave_Mechanical_in = sum(-(-P_out<=0).*(-P_out))*dt;
ave_Mechanical_out = sum((-P_out>=0).*(-P_out))*dt;
ave_inst_e = -(ave_Hydraulic_out-ave_Mechanical_out)/(ave_Hydraulic_in-ave_Mechanical_in)

perc_pos_power_in = ave_Hydraulic_in/(ave_Hydraulic_in+ave_Hydraulic_out)

