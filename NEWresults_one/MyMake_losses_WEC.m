RailFraction = 1;

load('Li_PI_one.mat')
V1 = myoutput.signals.values(:,11);
F1 = -myoutput.signals.values(:,17);
t = myoutput.time;
dt = mean(diff(t));
t_step = dt;
ACap1 = .2382;
ARod1 = ACap1;

% This code is used to change the drive cycle. Somethings may only be
% possible if power in changes sign (eg zero net flow).
% ind10 = find(t==100);
% V1(ind10:ind10+500) = max(V1);
% plot(t,V1)
% plot(t,cumsum(V1.*F1)*dt)


 
% This code was from Aditya. It was used to test whether this code was
% making the same decisions as adityas code was.
% load JCB_LongerDC.mat
% t_all = 1:length(Trenching.time);
% this_chunk = t_all;
% t = Trenching.time(this_chunk);
% t_step = mean(diff(t));
% F1 = -Trenching.boom_f(this_chunk);
% V1 = Trenching.boom_v(this_chunk);
% ACap1 = boom_Ac;
% ARod1 = boom_Ar;

% figure(50)
% subplot(2,1,1)
% plot(t,F1)
% legend('Force')
% subplot(2,1,2)
% plot(t,V1)
% legend('Velocity')
% 
% figure(51)
% plot(t,cumsum(F1.*V1)*dt)
% title('Power In')
% ylabel('Power (W)')
% xlabel('Time (s)')






%lambda1 = 0.25;

%% Main Pump Properties
%[Pmap_MP, Qmap_MP, Wmap_MP, DMP, CfMP, ChMP, CvMP, CsMP, CstMP] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');

xMP=1; % Frac Disp
wMP = 3000*2*pi/60; % Angular Velocity [rad/s]

%% HECM and MP properties

load AP_107cc_Sept16.mat

D=107e-6/(2*pi); %displacement [m^3/rad]

%% Load Drive Cycle/ Actuator Dimensions / HECM Dimensions

D_HECM1 = max(V1)*max([ACap1 ARod1])/(3000/60);   %assume 3000 rpm max speed
ScaleHECM1 = D_HECM1/(D*2*pi);
ScaleHECM1 = 45;
% HECM scale. calculated with ratio of max actuator flow and max map flow(?)

%PR = [0 1/2 1]; ScaleMaxT1 = 1; 
PR = [0, 1/6, 2/3, 1]; ScaleMaxT1 = 2.3;




nR = length(PR);
Frange1 = (ones(nR,1)*PR*ACap1-PR'*ones(1,nR)*ARod1); %forces that can be achieved by the pressure rails alone (at this stage in the code, PR = [0 1/3 1]
Frange1 = sort(unique(Frange1(:))); % order then in ascending order, and take out combos that repeat the forces

Pmax=max([max(-F1(-F1>0))/max(Frange1),min(-F1(-F1<0))/min(Frange1)]);
%Pmax=max([max(-F1>0)/max(Frange1),min(-F1<0)/min(Frange1)])

Pmax = 35e6;

if Pmax > 1,
    PR = PR*Pmax; 
    Frange1=Frange1*Pmax;
end
MaxT1_Act=max(diff(Frange1))/ARod1/2*ScaleHECM1*D*ScaleMaxT1; %e-motor torque limit T = deltaP*D

figure(1);
plot(t, -F1, t([1,end]),ones(2,1)*Frange1',t([1,end]),ones(2,1)*MaxT1_Act'); ylabel('Force (N)'); xlabel('Time (s)');grid

%% Given PR Selection and Directional Valve Indicator

DV_Ind = [1];

%% Find Pressure Differentials For HECMs, Pressures in the Cylinders (Depending on Which Combination of Rails/Valve and the Force)
[DV_Ind_HECM1, PRA1, PRB1] = ndgrid(DV_Ind,PR,PR); % Create a mesh grid of pressure to get all "combinations."

DV_Ind_HECM1 = DV_Ind_HECM1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA1 = PRA1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB1 = PRB1(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F1 = PRA1*ACap1-PRB1*ARod1;
MaxF1Diff = max(diff(sort(CPR_F1)))*RailFraction;% this line and the one above it dont seem to be used anywhere else

ones_t = ones(length(t),1);
ones_opt=ones(1,length(PRA1));
%% Find Angular Velocity of HECM Units @ Every Time


P1_Cyl = (ones_t*PRA1'*ACap1+F1*ones_opt)/ARod1;
DeltaP_HECM1 = P1_Cyl-ones_t*PRB1';

Q_HECM_1 = -repmat(V1,1,length(PRA1)).*ARod1; % velocity times ARod1 at every time (flow) - its the same for each PR combo

w1 = interp2(P_Map,ScaleHECM1*Q_Map,W_Map,DeltaP_HECM1,Q_HECM_1);
%w1 = interp2(P_diff,ScaleHECM1*Q_map,w_pump_map,DeltaP_HECM1,Q_HECM_1); %
%For minals map


% T1_Ideal = DeltaP_HECM1.*D*ScaleHECM1;
% T2_Ideal = DeltaP_HECM2.*D*ScaleHECM2;
% T3_Ideal = DeltaP_HECM3.*D*ScaleHECM3;
% T4_Ideal = DeltaP_HECM4.*D*ScaleHECM4;a


T1_Act = ScaleHECM1*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM1,w1);
%T1_Act = ScaleHECM1*interp2(P_diff,w_pump_map,T_map,DeltaP_HECM1,w1); % For minals map

%% Total HECM Hydraulic Power Loss
PLoss_Hydraulic_HECM = abs(T1_Act.*w1-DeltaP_HECM1.*Q_HECM_1);

%% Electrical Losses HECM 1


% Assume a constant electrical efficiency.
ElectricEff1 = .9; 

battery_power = (T1_Act.*w1>0).*(T1_Act.*w1)/ElectricEff1   +   (T1_Act.*w1<0).*(T1_Act.*w1)*ElectricEff1; % Energy leaving battery has a + sign

PLoss_Electric_HECM = (T1_Act.*w1>0).*abs(T1_Act.*w1*(1/ElectricEff1-1))...
                     +(T1_Act.*w1<0).*abs(T1_Act.*w1*(1-ElectricEff1)); 
                 


%% Find Main Pump + eGen Flows and Losses (Rail Losses)


wMP = 2000*2*pi/60; % Angular Velocity [rad/s]
T_PR = interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
Q_PR = interp2(P1_Mapping,w1rad_Mapping,Q1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
%T_PR = interp2(P_diff,w_pump_map,T_map,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP); % for minals map
%Q_PR = interp2(P_diff,w_pump_map,Q_map,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP); % for minals map

                                                             
MP_PumpEff = ElectricEff1*(Q_PR(1:(length(PR)-1)).*PR(2:end)')./(T_PR(1:(length(PR)-1))*wMP);
MP_MotorEff = ElectricEff1*T_PR(length(PR):end)*wMP./(-PR(2:end)'.*Q_PR(length(PR):end)); 

%MP_loss=(T_PR*wMP-Q_PR.*[PR(2:3)';-PR(2:3)'])./(Q_PR.*[PR(2:3)';PR(2:3)'])
R_P_loss = (1./MP_PumpEff-1);  % energy loss per unit volume of flow on the rail
R_M_loss = (1-MP_MotorEff);
%R_P_loss=[ 0.3361 ;0.2479; 0.2545];
%R_M_loss=[ 0.2795; 0.2042; 0.2059];

%R_P_loss=[ 0.2513; 0.2545];
%R_M_loss=[ 0.2085; 0.2059];

for k=1:length(PR)
    QR{k} = zeros(length(t),length(PRA1));
    ind = find(PRA1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+V1*ACap1*ones(1,length(ind));
    ind = find(PRB1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+Q_HECM_1(:,ind);
end

%% Sum All Losses/Costs ***************************************************

% Torque limit
PLoss_Electric_HECM(abs(T1_Act) > MaxT1_Act) = inf;

% Cavitation limit
PLoss_Hydraulic_HECM(P1_Cyl < -1e5) = inf;

HECMLosses = PLoss_Hydraulic_HECM +PLoss_Electric_HECM;
%battery_power; 
%%
figure(), plot(t,HECMLosses), xlim([140 160])

