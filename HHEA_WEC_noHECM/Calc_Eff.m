% Calculate the efficiency of the main pump
% This code is called in WEC_DP_v7

% It takes as inputs:
    % Factional displacements
    % Pressure rails
    % Time at each rail
    % size of motor

% This model for losses comes from a paper written by someone else % Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
% e-map is just constant 90% efficiency

fracDisp = [chi_HT chi_MT chi_HL];
%fracDisp = linspace(0,1,20);
%fracDisp = 1;
%wMP = linspace(0,315,20);
%wMP = 16;
deltaP = [PR(4) PR(3) (PR(4)-PR(2))];
%deltaP = linspace(0,PR(4),20);
Scale = D*2*pi*1e6/107; % How much larger is this pump than the 107 cc pump for which this model was created?
[a,b] = meshgrid(fracDisp,deltaP);
%[a,b] = meshgrid(wMP,deltaP);

for i = 1:length(b) % loop through deltaP
    for j = 1:length(a) %loop through fracDisp or wMP
        [P_in(i,j), P_out(i,j)] = motor_map_func(deltaP(i),fracDisp(j),wMP,Scale);
        %[P_in(i,j), P_out(i,j)] = motor_map_func(deltaP(i),fracDisp,wMP(j),Scale);
    end
end


e = P_out ./ P_in ;
%figure, contour(a,b,e,.1:.1:1,'ShowText','on')

total_MP_e = ( t_HT*e(1,1) + t_MT*e(2,2) + t_HL*e(3,3)) / T
return
%%
%[in, out] = motor_map_func(35e6,fracDisp,314,Scale);
%out/in

function [P_in P_out] = motor_map_func(deltaP,fracDisp,wMP,Scale)
% This model for losses comes from a paper written by someone else % Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
%% Fluid Properties
mu=(32e-6)*870;
B = 1.7e9;
rho = 870;
%% Manufacturer 107cc/rev
% Variable Displacement Axial Piston, 107 cc/rev (Pourmovahed et al. 1992b)
d = (107*100^-3)/(2*pi); % m^3/rad
% Torques Loss Constants
Cf =  53.7e-3;
Ch = 53.6;
Cv = 23.5e3;
% Flow Loss Constants
Cs = 4.26e-9;
Cst = 0*1e-5;

Q_Ideal = fracDisp*d*wMP*Scale;
QLoss = Scale*abs(d*Cs*(deltaP)/mu) + Scale*abs(fracDisp*d*wMP*(deltaP)/B) + Scale*abs(d^(2/3)*Cst*(2*(deltaP)/rho)^.5);
Q_Act = Q_Ideal + QLoss; % only valid for motoring

T_Ideal = deltaP*d*fracDisp*Scale;
TLoss = Scale*(  abs(d*Cv*mu*wMP) + abs(d*(deltaP)*Cf) + abs(fracDisp*Ch*wMP^2*rho*d^(5/3)/2)  );
T_Act = T_Ideal - TLoss; % only valid for motoring


% Power in
P_in = deltaP*Q_Act;

% Power out with 90% effiency
P_out = .9*wMP*T_Act;
end