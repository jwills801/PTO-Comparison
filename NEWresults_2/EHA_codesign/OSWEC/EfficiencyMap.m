%% Pump Motor Efficiency Mapping

%% Fluid Properties

mu=(32*10^-6)*870;
B = 1.7*10^9;
rho = 870;

%% Define Pump Constants

%%

%% UW Madison
% % Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
%     D = 107; % cc/rev
%     d = (D*100^-3)/(2*pi); % m^3/rad
% % Torques Loss Constants
%     Cf = 4.8*10^-3;
%     Ch = 0;
%     Cv = 153*10^3;
% % Flow Loss Constants
%     Cs = 1.04*10^-9;
%     Cst = 1.20*10^-5;
    
%% Manufacturer 107cc/rev
% Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
    D = 107; % cc/rev
    d = (D*100^-3)/(2*pi); % m^3/rad
% Torques Loss Constants
    Cf = 53.7*10^-3;
    Ch = 53.6;
    Cv = 23.5*10^3;
% Flow Loss Constants
    Cs = 4.26*10^-9;
    Cst = 0*10^-5;

%% Define Pump Parameters to Study
% (fractional displacement/angular velocity/pressure differential)

% Angular Velocity
Wrpm = 3000;
w = Wrpm.*(2*pi/60); % radians per second

% Fractional Displacement Range
fracDisp = linspace(-1,1,50);

% Pressure Differential (Pressure actuator - Pressure Rail)
deltaP = -40e6:1e6:40e6;
% deltaP = [.5 1]*35e6;
%% Find Torque and Flow Losses

QLoss = zeros(length(fracDisp),length(deltaP));
Q_Ideal = zeros(length(fracDisp),length(deltaP));
Q_Act = zeros(length(fracDisp),length(deltaP));

TLoss = zeros(length(fracDisp),length(deltaP));
T_Ideal = zeros(length(fracDisp),length(deltaP));
T_Act = zeros(length(fracDisp),length(deltaP));
Q_Efficiency = zeros(length(fracDisp),length(deltaP));
T_Efficiency = zeros(length(fracDisp),length(deltaP));

Quad = zeros(length(fracDisp),length(deltaP));


for i = 1:length(fracDisp)
    for j = 1:length(deltaP)
         
%         Q_Ideal(i,j) = w(i)*d*fracDisp;
%         T_Ideal(i,j) = deltaP(j)*d*fracDisp;
%         QLoss(i,j) = abs(d*Cs*(deltaP(j))/mu) + abs(fracDisp*d*w(i)*(deltaP(j))/B) + abs(d^(2/3)*Cst*(2*(deltaP(j))/rho)^.5);
%         TLoss(i,j) = abs(d*Cv*mu*w(i)) + abs(d*(deltaP(j))*Cf) + abs(fracDisp*Ch*w(i)^2*rho*d^(5/3)/2);
        
        Q_Ideal(i,j) = fracDisp(i)*d*w*Scale;
        T_Ideal(i,j) = deltaP(j)*d*fracDisp(i)*Scale;
        QLoss(i,j) = Scale*abs(d*Cs*(deltaP(j))/mu) + Scale*abs(fracDisp(i)*d*w*(deltaP(j))/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP(j))/rho)^.5);
        TLoss(i,j) = Scale*(  abs(d*Cv*mu*w) + abs(d*(deltaP(j))*Cf) + abs(fracDisp(i)*Ch*w^2*rho*d^(5/3)/2)  );
%         Q_Act(i,j) = Q_Ideal-sign(w(i)*deltaP(j))*QLoss(i,j)
%         T_Act(i,j) = T_Ideal+sign(w(i)*deltaP(j))*TLoss(i,j)
        
        if deltaP(j)>=0 && w >=0
            Quad(i,j) = 1;
            Q_Act(i,j) = Q_Ideal(i,j) - QLoss(i,j);
            T_Act(i,j) = T_Ideal(i,j) + TLoss(i,j);
        
        elseif deltaP(j)<0 && w >=0
            Quad(i,j) = 4;
            Q_Act(i,j) = Q_Ideal(i,j) + QLoss(i,j);
            T_Act(i,j) = T_Ideal(i,j) + TLoss(i,j);
        
        elseif deltaP(j)>=0 && w <0
            Quad(i,j) = 2;
            Q_Act(i,j) = Q_Ideal(i,j) - QLoss(i,j);
            T_Act(i,j) = T_Ideal(i,j) - TLoss(i,j);
            
        else
            Quad(i,j) = 3;
            Q_Act(i,j) = Q_Ideal(i,j) + QLoss(i,j);            
            T_Act(i,j) = T_Ideal(i,j) - TLoss(i,j);
        end
        
        if abs(Q_Act(i,j)) >= abs(Q_Ideal(i,j))
            Q_Efficiency(i,j) = Q_Ideal(i,j)/Q_Act(i,j);
        else
            Q_Efficiency(i,j) = Q_Act(i,j)/Q_Ideal(i,j);
        end
        
       if abs(T_Act(i,j)) >= abs(T_Ideal(i,j))
           T_Efficiency(i,j) = T_Ideal(i,j)/T_Act(i,j);
       else
           T_Efficiency(i,j) = T_Act(i,j)/T_Ideal(i,j);
       end
        
     end
end
% QEfficiency = Q_Act./Q_Ideal;
Q_Efficiency = Q_Efficiency';
T_Efficiency = T_Efficiency';
Total_Efficiency = (Q_Efficiency.*T_Efficiency);

%% Plot Contours
%levels = 0:.05:.95;



% % Volumetric
% figure
% %levels = [.995 .99 .985 .98 .975 .97];
% levels = 0:.1:1;
% contour(w,(deltaP),Q_Efficiency,levels,'ShowText','on')
% grid on,
% ylabel('Delta Pressure [MPa]')
% xlabel('Angular Velocity [rad/s]')
% %set(gcf, 'Position', get(0, 'Screensize'));
% 
% % Torque
% figure
% %levels = .9:.01:.99;
% contour(w,(deltaP),T_Efficiency,levels,'ShowText','on')
% grid on
% ylabel('Delta Pressure [MPa]')
% xlabel('Angular Velocity [rad/s]')
% %set(gcf, 'Position', get(0, 'Screensize'));

% Overall
levels = .1:.05:1.05;
levels = .25:.25:1;
contour(fracDisp,(deltaP),Total_Efficiency,levels,'ShowText','on')
grid on
ylabel('$\Delta$ P [MPa]','interpreter','latex')
xlabel('Fractional Displacement [rad/s]','interpreter','latex')
% set(gcf, 'Position', get(0, 'Screensize'));