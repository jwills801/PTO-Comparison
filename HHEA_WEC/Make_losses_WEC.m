tic
RailFraction = 1;

load('wec-PI_out_A_alternate.mat')
F1 = -output.controller.force;
V1 = output.controller.velocity;
t = output.controller.time;
dt = mean(diff(t));

%lambda1 = 0.25;

%% Main Pump Properties
%[Pmap_MP, Qmap_MP, Wmap_MP, DMP, CfMP, ChMP, CvMP, CsMP, CstMP] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');

xMP=1; % Frac Disp
wMP = 2000*2*pi/60; % Angular Velocity [rad/s]

%% HECM and MP properties
disp('Pump maps')
load AP_107cc_Sept16.mat

D=107e-6/(2*pi); %displacement [m^3/rad]

%% Load Drive Cycle/ Actuator Dimensions / HECM Dimensions
disp('Drive cycle')
% load CNH Wheel Loader Drive Cycle
% *** NOTE that the forces are flipped. This is because the data file
% provided forces with a positive sign when acting in the direction from
% rod to cap end. The following code assumed a positive force when acting
% from cap to rod end on the cylinder.
% Piston areas
ACap1 = 0.0378; % This used to be 0.0378 
ARod1 = 0.0378;

disp('HECM pump sizes in cc/rev')
D_HECM1 = max(V1)*max([ACap1 ARod1])/(3000/60);   %assume 3000 rpm max speed
ScaleHECM1 = D_HECM1/(D*2*pi);

% HECM scale. calculated with ratio of max actuator flow and max map flow(?)
ScaleMaxT1 = 80; % 4.61

PR = [0 1/3 1]; % Define an array of pressures
%PR = [0, 1/6, 2/3, 1];

nR = length(PR);
Frange1 = (ones(nR,1)*PR*ACap1-PR'*ones(1,nR)*ARod1); %forces that can be achieved by the pressure rails alone (at this stage in the code, PR = [0 1/3 1]
Frange1 = sort(unique(Frange1(:))); % order then in ascending order, and take out combos that repeat the forces

Pmax=max([max(-F1(-F1>0))/max(Frange1),min(-F1(-F1<0))/min(Frange1)]);
%Pmax=max([max(-F1>0)/max(Frange1),min(-F1<0)/min(Frange1)])

Pmax = 35e6;
disp('Updated rail pressures')
if Pmax > 1,
    PR = PR*Pmax; 
    Frange1=Frange1*Pmax;
end
MaxT1_Act=max(diff(Frange1))/ARod1/2*ScaleHECM1*D*ScaleMaxT1; %e-motor torque limit T = deltaP*D

figure(1);
plot(t, -F1, t([1,end]),ones(2,1)*Frange1',t([1,end]),ones(2,1)*MaxT1_Act'); ylabel('Force1 N'); xlabel('Time - s');grid

%% Given PR Selection and Directional Valve Indicator

DV_Ind = [1];

%% Find Pressure Differentials For HECMs, Pressures in the Cylinders (Depending on Which Combination of Rails/Valve and the Force)
disp('Setup')
[DV_Ind_HECM1, PRA1, PRB1] = ndgrid(DV_Ind,PR,PR); % Create a mesh grid of pressure to get all "combinations."

DV_Ind_HECM1 = DV_Ind_HECM1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA1 = PRA1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB1 = PRB1(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F1 = PRA1*ACap1-PRB1*ARod1;
MaxF1Diff = max(diff(sort(CPR_F1)))*RailFraction;% this line and the one above it dont seem to be used anywhere else

ones_t = ones(length(t),1);
ones_opt=ones(1,length(PRA1));
%% Find Angular Velocity of HECM Units @ Every Time
disp('HECM Hydraulics')

disp('HECM Hydraulics - P')

P1_Cyl = (ones_t*PRA1'*ACap1+F1*ones_opt)/ARod1;
DeltaP_HECM1 = P1_Cyl-ones_t*PRB1';

disp('HECM Hydraulics - Q')

Q_HECM_1 = -repmat(V1,1,length(PRA1)).*ARod1; % velocity times ARod1 at every time (flow) - its the same for each PR combo

disp('HECM Hydraulics - w')

w1 = interp2(P_Map,ScaleHECM1*Q_Map,W_Map,DeltaP_HECM1,Q_HECM_1);

% T1_Ideal = DeltaP_HECM1.*D*ScaleHECM1;
% T2_Ideal = DeltaP_HECM2.*D*ScaleHECM2;
% T3_Ideal = DeltaP_HECM3.*D*ScaleHECM3;
% T4_Ideal = DeltaP_HECM4.*D*ScaleHECM4;
disp('HECM Hydraulics - T')

T1_Act = ScaleHECM1*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM1,w1);

%% Total HECM Hydraulic Power Loss
disp('HECM Hydraulics - Loss')

PLoss_Hydraulic_HECM = abs(T1_Act.*w1-DeltaP_HECM1.*Q_HECM_1);

%% Electrical Losses HECM 1
disp('HECM electric')

% Assume a constant electrical efficiency.
ElectricEff1 = 0.9; 

battery_power = (T1_Act.*w1>0).*(T1_Act.*w1)/ElectricEff1   +   (T1_Act.*w1<0).*(T1_Act.*w1)*ElectricEff1; % Energy leaving battery has a + sign

PLoss_Electric_HECM = (T1_Act.*w1>0).*abs(T1_Act.*w1*(1/ElectricEff1-1))...
                     +(T1_Act.*w1<0).*abs(T1_Act.*w1*(1-ElectricEff1)); 
                 
%% Find Main Pump + eGen Flows and Losses (Rail Losses)
disp('Main pump losses');

wMP = 2000*2*pi/60; % Angular Velocity [rad/s]
T_PR = interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
Q_PR = interp2(P1_Mapping,w1rad_Mapping,Q1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
                                                             
MP_PumpEff = ElectricEff1*(Q_PR(1:(length(PR)-1)).*PR(2:end)')./(T_PR(1:(length(PR)-1))*wMP);
MP_MotorEff = ElectricEff1*T_PR(length(PR):end)*wMP./(-PR(2:end)'.*Q_PR(length(PR):end)); 

%MP_loss=(T_PR*wMP-Q_PR.*[PR(2:3)';-PR(2:3)'])./(Q_PR.*[PR(2:3)';PR(2:3)'])
R_P_loss = (1./MP_PumpEff-1);  % energy loss per power
R_M_loss = (1-MP_MotorEff);
%R_P_loss=0*R_P_loss; %debug
%R_M_loss=0*R_M_loss;

for k=1:length(PR)
    QR{k} = zeros(length(t),length(PRA1));
    ind = find(PRA1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+V1*ACap1*ones(1,length(ind));
    ind = find(PRB1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+Q_HECM_1(:,ind);
end

%% Sum All Losses/Costs ***************************************************
disp('Constraints')
% Torque limit
PLoss_Electric_HECM(abs(T1_Act) > MaxT1_Act) = inf;
% Cavitation limit
%PLoss_Hydraulic_HECM(P1_Cyl < -1e5) = inf;

HECMLosses = PLoss_Hydraulic_HECM +PLoss_Electric_HECM;
%battery_power; 
%%
toc
