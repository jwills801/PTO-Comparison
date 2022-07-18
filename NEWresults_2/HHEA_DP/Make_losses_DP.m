tic
% The HECM is on the rod side of each actuator.

ACap1 = .2382; ARod1 = ACap1;

% HECM and MP properties
load ../AP_107cc_Sept16.mat; % Pump map

D_HECM1 = max(abs(V1))*ARod1/(maxRPM/60);

ScaleHECM1 = D_HECM1/(Disp*2*pi);

nR = length(PR);
Frange1 = (ones(nR,1)*PR*ACap1-PR'*ones(1,nR)*ARod1);
Frange1 = sort(unique(Frange1(:)));


figure, plot(t, -F1, t([1,end]),ones(2,1)*Frange1'); ylabel('Force (N)'); xlabel('Time (s)'); grid


%% Find Pressure Differentials For HECMs, Pressures in the Cylinders (Depending on Which Combination of Rails/Valve and the Force)
[PRA1,PRB1] = ndgrid(PR,PR); % Create a mesh grid of pressure to get all "combinations."

PRA1 = PRA1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB1 = PRB1(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F1 = PRA1*ACap1-PRB1*ARod1;

ones_t = ones(length(t),1);
ones_opt=ones(1,length(PRA1));
%% Find Angular Velocity of HECM Units @ Every Time
P1_Cyl = (ones_t*PRA1'*ACap1+F1*ones_opt)/ARod1; % pressure in rod side
DeltaP_HECM1 = P1_Cyl-ones_t*PRB1';

Q_HECM_1 = -V1*ARod1;

w1 = interp2(P_Map,ScaleHECM1*Q_Map,W_Map,DeltaP_HECM1,Q_HECM_1*ones_opt);
T1_Act = ScaleHECM1*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM1,w1);

PLoss_Hydraulic_HECM = T1_Act.*w1-DeltaP_HECM1.*(Q_HECM_1*ones_opt);

% Torque limit
[MaxT1_Act] = min_torque_limit(T1_Act,P1_Cyl,t_c,spacing);


%% Electrical Losses HECM

switch Emap
    case 'Mani' 
        load ../E_Map_ManikantaPalantla_3point2kW.mat
        sfw = 1.15;
        speedscale = sfw*maxRPM*(2*pi/60)/max(abs(w_m(:)));
        Tscale1 = MaxT1_Act/max(abs(T_m(:)));
        
        PLoss_Electric_HECM = interp2(speedscale*w_m,Tscale1*T_m,speedscale*Tscale1*ELoss_m,w1,T1_Act);

    case 'constantEff'
        ElectricEff1 = 0.9;

        PLoss_Electric_HECM = (T1_Act.*w1>0).*abs(T1_Act.*w1*(1/ElectricEff1-1))...
                     +(T1_Act.*w1<0).*abs(T1_Act.*w1*(1-ElectricEff1));
                 
    otherwise
        disp('Please define Electric Motor Map');
        pause;
end

battery_power = (T1_Act.*w1) + PLoss_Electric_HECM;
               
%% Find Main Pump Flows and Losses (Rail Losses)
wMP = 2000*2*pi/60; % Angular Velocity [rad/s]
T_PR = interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
Q_PR = interp2(P1_Mapping,w1rad_Mapping,Q1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);

MP_PumpEff = (Q_PR(1:(length(PR)-1)).*PR(2:end)')./(T_PR(1:(length(PR)-1))*wMP);
MP_MotorEff = T_PR(length(PR):end)*wMP./(-PR(2:end)'.*Q_PR(length(PR):end));

% Inlcuding a constant electrical efficiency of 90%
R_P_loss = (1/1./MP_PumpEff-1);  % loss per power
R_M_loss = (1-1*MP_MotorEff);

%% 
for k=1:length(PR)
    QR{k} = zeros(length(t),length(PRA1));
    ind = find(PRA1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+V1*ACap1*ones(1,length(ind));
    ind = find(PRB1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+Q_HECM_1*ones(1,length(ind));
end

%% Constraints
PLoss_Electric_HECM(abs(T1_Act) >= MaxT1_Act) = inf;
% Cavitation limit
PLoss_Hydraulic_HECM(P1_Cyl < -1e5) = inf;
HECMLosses = PLoss_Hydraulic_HECM+PLoss_Electric_HECM;

%%
disp(['Making Elec and Hydr losses took ' num2str(toc) ' seconds'])