function [TotalCost, M, I, genP, choices, ScaleHECM1, params] = HHEA_pto(F,v,t)
tic
% HHEA version from 7/30/2019

startingInd = 1;
spacing = 1;
% lambda1 = 0.244643402099609; % lambda for JCB cycle, One Pressure Rail + Trenching
lambda1 = 0.25; % using equivalent lambda for WEC-Sim case


%% Fluid Properties
mu = (32e-6)*870; % Dynamic Viscosity []
B = 1.7e9; % Bulk Modulus []
rho = 870; % Density [kg/m^3]


%% Main Pump Properties
[Pmap_MP, Qmap_MP, Wmap_MP, DMP, CfMP, ChMP, CvMP, CsMP, CstMP] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');
xMP = 1; % Main Pump Fracional Displacement
wMP = 2000*2*pi/60; % Angular Velocity [rad/s]


%% HECM Pump 1 Properties
[Pmap1, Qmap1, Wmap1, D1, Cf1, Ch1, Cv1, Cs1, Cst1] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');
x1 = 1; % HECM Fractional Displacement

%%
load('AP_107cc_M_FullDisp_TESTpumpINFOfunction.mat')

%% Load Drive Cycle/ Actuator Dimensions / HECM Dimensions

% Drive Cycle *****************************************************
% Actuator load force, defined as positive when trying to force extension
% Actuator velocity, defined as positive in extension
%  _____________________________
% |           |
% |           |----------------------------   -> +Force, F1
% |___________|_________________              -> +Velocity, V1
%          Actuator1

t_step = t(2)-t(1); % equal to 0.1 for WEC-Sim cases
F1 = F; % actuator force for HECM 1
V1 = v; % actuator velocity for HECM 1

% Piston areas
ACap1 = 0.0378;
ARod1 = 0.0378;

% HECM scale. calculated with ratio of max actuator flow and max map flow
ScaleHECM1 = max(abs(v))*ACap1/2000*60*1e6 / 100;


%% Given PR Selection and Directional Valve Indicator
PR = [0 1/2 1]*35e6; % Pressure rail pressures (tank, low, high)
DV_Ind = [1 0]; % index to determine if HECM connected to cap or rod side


%% Find Pressure Differentials For HECMs, Pressures in the Cylinders (Depending on Which Combination of Rails/Valve and the Force)
[DV_Ind_HECM1, PRA1, PRB1] = ndgrid(DV_Ind,PR,PR); % Create a mesh grid of pressure to get all "combinations."

DV_Ind_HECM1 = DV_Ind_HECM1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA1 = PRA1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB1 = PRB1(:); % Spread the combinations of pressure rails and directional valves into a vector.

PA1 = PRA1; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on B Side)
PB1 = PRB1; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on A Side)

CPR_A1 = repmat(PA1',length(t),1);
CPR_B1 = repmat(PB1',length(t),1);

P1_Cyl = zeros(length(t),length(PRA1)); % Pressure on actuator side connected to HECM (regardless if cap/rod)
DeltaP_HECM1 = zeros(length(t),length(PRA1));

for i = 1:size(PRA1,1)
    if DV_Ind_HECM1(i) == 0 
        P1_Cyl(:,i) = (PB1(i)*ARod1-F1)/ACap1;
        DeltaP_HECM1(:,i) = P1_Cyl(:,i)-PRA1(i);
    else
        P1_Cyl(:,i) = (PA1(i)*ACap1+F1)/ARod1;
        DeltaP_HECM1(:,i) = P1_Cyl(:,i)-PRB1(i);
    end
end


%% Find Angular Velocity of HECM Units @ Every Time
V1_Array = repmat(V1,1,length(PRA1)); % Repeat the desired velocity matrix for next calculation
Q1_A = V1_Array.*ACap1; % The flow taken from rail connected to the A side.(HECM1)
Q1_B = -V1_Array.*ARod1; % The flow taken from rail connected to the B side.(HECM1)

Q_HECM_1 = zeros(length(t),length(PRA1));

for i = 1:length(DV_Ind_HECM1)
    if DV_Ind_HECM1(i) == 0
        Q_HECM_1(:,i) = Q1_A(:,i);
    else
        Q_HECM_1(:,i) = Q1_B(:,i);
    end

end

w1 = interp2(Pmap1,ScaleHECM1*Qmap1,Wmap1,DeltaP_HECM1,Q_HECM_1);
% T1_Ideal = DeltaP_HECM1.*D1*ScaleHECM1;
T1_Act = ScaleHECM1*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM1,w1);


%% Total HECM Hydraulic Power Loss
PLoss_Hydraulic_HECM1_check = abs(T1_Act.*w1-DeltaP_HECM1.*Q_HECM_1);
PLoss_Hydraulic_HECM = PLoss_Hydraulic_HECM1_check;


%% Electrical Losses HECM 1
% Assume a constant electrical efficiency.
ElectricEff1 = .90;

w1_Generator = w1; % Grabs the angular velocity of HECM
w1_Generator(sign(T1_Act.*w1)>=0) = 0; % Throws away values when Motoring
w1_Motor = w1; % Grabs the angular velocity of HECM
w1_Motor(sign(T1_Act.*w1)<0) = 0; % Throws away values when Generating

PLoss1_Electric_Generator = abs(T1_Act.*w1_Generator*(1-ElectricEff1));
PLoss1_Electric_Motor = abs(T1_Act.*w1_Motor*(1/ElectricEff1-1));

PLoss1_Electric = PLoss1_Electric_Generator+PLoss1_Electric_Motor;

% Electrical power /output/ of the HECM (elec power to shore, + if gen, - if mot)
powerHECM = -1*T1_Act.*w1_Generator - PLoss1_Electric_Generator + ...
    -1*T1_Act.*w1_Motor + PLoss1_Electric_Motor;


%% Total Electric Power Loss from HECMs
PLoss_Electric_HECM = PLoss1_Electric;


%% Find Main Pump & Main Generator Losses (Rail Losses)
MainPumpEfficiency = .885;
totalEff = MainPumpEfficiency * ElectricEff1; % 0.7965 with 90% electric & 88.5% MP

% rail flows
QR1 = Q1_A.*(CPR_A1==PR(1))+Q1_B.*(CPR_B1==PR(1));
QR2 = Q1_A.*(CPR_A1==PR(2))+Q1_B.*(CPR_B1==PR(2));
QR3 = Q1_A.*(CPR_A1==PR(3))+Q1_B.*(CPR_B1==PR(3));

% main pump losses
MP_PumpEff = totalEff;
MP_MotorEff = totalEff;
R_P_loss = (1./MP_PumpEff-1);  % loss per power for pumping
R_M_loss = (1-MP_MotorEff);   % loss per power for motoring
MPLosses = abs(QR2)*PR(2)*R_M_loss + abs(QR3)*PR(3)*R_M_loss; % use abs(Q) before optimization, or actual flow after optimizing

powerMP = totalEff*-1*(Q1_A.*CPR_A1+Q1_B.*CPR_B1); % total output electrical power of the main pump


%% Sum All Losses/Costs ***************************************************
TotalLosses = PLossMP+PLoss_Hydraulic_HECM+PLoss_Electric_HECM;
Cost = PLossMP+PLoss_Hydraulic_HECM+PLoss_Electric_HECM;

% Set losses for combinations that produce cavitation to inf so that they
% are not picked!
TotalLosses(P1_Cyl<0) = inf;
Cost(P1_Cyl<0) = inf;


%% Find Minimum
% Find the combination with the least power loss. "I" represents which
% combination produces the minimum power loss and M is the power loss for
% that combination.
[M, I] = min(Cost,[],2);
TotalCost = sum(M)*t_step*spacing;

genP1 = F1.*V1 - M;
genP2 = zeros(length(t),1);
genP3 = zeros(length(t),1);
phecm = zeros(length(t),1);
pmp = zeros(length(t),1);
qr1 = zeros(length(t),1);
qr2 = zeros(length(t),1);
qr3 = zeros(length(t),1);
choices = zeros(length(t),3);

for i=1:length(t)
      % actual HECM and Mp electrical power output
      powerHECM(i) = powerHECM(i,I(i));
      pmp(i) = powerMP(i,I(i));
      
      % actual choices of directional valve and pressure rails
      choices(i,1) = DV_Ind_HECM1(I(i));
      choices(i,2) = PRA1(I(i));
      choices(i,3) = PRB1(I(i));
      
      % actual hecm and MP flow rates
      hecmFlow = Q_HECM_1(i,I(i));
      mpFlow = hecmFlow;
      
      % actual flow rates for each rail
      qr1(i) = QR1(i,I(i));
      qr2(i) = QR2(i,I(i));
      qr3(i) = QR3(i,I(i));
end

% % can also use this statement to calculate MP losses after optimizing
% % check sign, +flow = motoring or pumping?
% mp_output = (qr1(qr1>0)*pr1 + qr2(qr2>0)*pr2)*totalEff + ...
%     (qr1(qr1<0)*pr1 + qr2(qr2<0)*pr2)*(1/totalEff);

genP = phecm + pmp; % total generated electrical power output

mpDisplacement = inf; % what should this be for this scale?
mpAngVel = mpFlow / mpDisplacement;
mpTorque = pmp ./ mpAngVel;

toc

%% save outputs to a struct
params = struct();
params.actForce = F1; % actuator force
params.actVel = V1; % actuator velocity
params.gen_power = genP; % total electrical power generated
params.mp_power = pmp; % main pump electrical power output
params.hecm_power = phecm; % HECM electrical power output
params.flow = hecmFlow; % flow rate through HECM 
params.mpTorque = mpTorque; % main pump torque -- needs correct displacement
params.mpAngVel = mpAngVel; % main pump angular velocity -- needs correct displacement


%% Optional - plot power loss components
if(false)
mp = zeros(1,length(t));
hyd = zeros(1,length(t));
elec = zeros(1,length(t));
gen = zeros(1,length(t));
qp = zeros(1,length(t));
fv = zeros(1,length(t));
for i=1:length(t)
    hyd(i) = PLoss_Hydraulic_HECM(i,I(i));
    elec(i) = PLoss_Electric_HECM(i,I(i));
    mp(i) = PLossMP(i,I(i));
    gen(i) = PLossMG;
    fv(i) = F1(i)*V1(i);
end

figure()
plot(t,fv,t,genP,t,M,t,mp,t,hyd,t,elec);
legend('Actuator intput power', 'total electrical power output',...
    'Total power loss', 'main pump power loss', 'HECM hydraulic power loss',...
    'HECM electrical power loss');
xlabel('Time (s)');
ylabel('Power losses (W)');
title('Power loss components');
end

