%% For Electrified HHEA

% tic
RailFraction = 1;
startingInd = 1;
spacing = 1;

% ScaleMaxT = .6;
% ScaleMaxT1 = ScaleMaxT;
% ScaleMaxT2 = ScaleMaxT;

%PR = [0 0.2 0.8 1]; % Define an array of pressures

%% HECM and MP properties
disp('Pump maps')
%load AP_107cc_Sept16.mat

D=Disp_28cc; %displacement [m^3/rad]  %addressing from the 28point1.... file

find_size;
D_HECM4 = 18.1*(1/100)^3    %rotary actuator sizing

disp('HECM pump sizes in cc/rev')
[D_HECM1,D_HECM2, D_HECM3, D_HECM4]*1e6
ScaleHECM1 = D_HECM1/(D*2*pi);
ScaleHECM2 = D_HECM2/(D*2*pi);
ScaleHECM3 = D_HECM3/(D*2*pi);
ScaleHECM4 = D_HECM4/(D*2*pi);
%% Load Drive Cycle/ Actuator Dimensions / HECM Dimensions
disp('Drive cycle')
% load CNH Wheel Loader Drive Cycle
% *** NOTE that the forces are flipped. This is because the data file
% provided forces with a positive sign when acting in the direction from
% rod to cap end. The following code assumed a positive force when acting
% from cap to rod end on the cylinder.

% Key:
% LinearActuator 1: Boom
% LinearActuator 2: Stick
% LinearActuator 3: Bucket
% RotaryActuator 4: Swing



%% Jackson Commented the inended stuff out
% % Drive Cycle***************************************************
% cycle = LC;
% cycle = SC;
                %t_all = 1:spacing:length(cycle.time);
% t_all = 1:spacing:100;
%                this_chunk = t_all;

%               t = cycle.time(this_chunk);
                %t_step = mean(diff(t));
                 %F1 = -cycle.boom_f(this_chunk);
               %V1 = cycle.boom_v(this_chunk);
            %ACap1 = boom_Ac;
             %ARod1 = boom_Ar;

%% Jacksons Changes
load wec-PI_out_A_alternate.mat
F1 = -output.controller.force;
V1 = output.controller.velocity;
t = output.controller.time;
t_step = mean(diff(t));

t_all = 1:1:length(t);
this_chunk = t_all;



ACap1 = 0.0578; % This used to be 0.0378 
ARod1 = 0.0378;



F2 = F1;
V2 = V1;
ACap2 = ACap1;
ARod2 = ARod1;

F3 = F1;
V3 = V1;
ACap3 = ACap1;
ARod3 = ARod1;


%%

%F2 = -cycle.arm_f(this_chunk);
%V2 = cycle.arm_v(this_chunk);
%ACap2 = arm_Ac;
%ARod2 = arm_Ar;

%F3 = -cycle.bkt_f(this_chunk);
%V3 = cycle.bkt_v(this_chunk);
%ACap3 = bkt_Ac;
%ARod3 = bkt_Ar;

T4 = -cycle.sw_t(this_chunk)/Swing_ratio;
W4 = cycle.sw_w(this_chunk)*Swing_ratio;

nR = length(PR);
Frange1 = (ones(nR,1)*PR*ACap1-PR'*ones(1,nR)*ARod1);
Frange1 = sort(unique(Frange1(:)));
Frange2 = (ones(nR,1)*PR*ACap2-PR'*ones(1,nR)*ARod2);
Frange2 = sort(unique(Frange2(:)));
Frange3 = (ones(nR,1)*PR*ACap3-PR'*ones(1,nR)*ARod3);
Frange3 = sort(unique(Frange3(:)));
Trange4 = sort(unique(ones(nR,1)*PR-PR'*ones(1,nR)))*ScaleHECM4*D;

% Pmax=max([max(-F1(-F1>0))/max(Frange1),min(-F1(-F1<0))/min(Frange1),...
%     max(-F2(-F2>0))/max(Frange2),min(-F2(-F2<0))/min(Frange2)])
%Pmax=max([max(-F1>0)/max(Frange1),min(-F1<0)/min(Frange1),...
%    max(-F2>0)/max(Frange2),min(-F2<0)/min(Frange2)])



% Pmax = 29e6;
Pmax = 1; %for LRP = 0.09*Pmax (based on 'finding_Pmax' code)
disp('Updated rail pressures')
if Pmax > 1,
    PR = PR*Pmax; 
    Frange1=Frange1*Pmax;
    Frange2=Frange2*Pmax;
    Frange3=Frange3*Pmax;
    Trange4=Trange4*Pmax;
end; PR

% MaxT1_Act = 101;
% MaxT2_Act = 74;
% MaxT1_Act=max(diff(Frange1))/ARod1/2*ScaleHECM1*D*ScaleMaxT1;
% MaxT2_Act=max(diff(Frange2))/ARod2/2*ScaleHECM2*D*ScaleMaxT2;
% MaxT1_Act=max(diff(Frange1))/ARod1*ScaleHECM1*D*ScaleMaxT1;
% MaxT2_Act=max(diff(Frange2))/ARod2*ScaleHECM2*D*ScaleMaxT2;

% figure(1);
% subplot(221); plot(t, -F1, t([1,end]),ones(2,1)*Frange1'); ylabel('Force1 N'); xlabel('Time - s');grid
% subplot(222); plot(t, -F2, t([1,end]),ones(2,1)*Frange2'); ylabel('Force2 N'); xlabel('Time - s');grid
% subplot(223); plot(t, -F3, t([1,end]),ones(2,1)*Frange3'); ylabel('Force3 N'); xlabel('Time - s');grid
% subplot(224); plot(t, -T4, t([1,end]),ones(2,1)*Trange4'); ylabel('Torque4 Nm'); xlabel('Time - s');grid

%% Given PR Selection and Directional Valve Indicator

DV_Ind = [1];

%% Find Pressure Differentials For HECMs, Pressures in the Cylinders (Depending on Which Combination of Rails/Valve and the Force)
disp('Setup')
[DV_Ind_HECM1, PRA1, PRB1, DV_Ind_HECM2, PRA2, PRB2, DV_Ind_HECM3, PRA3, PRB3, PRA4, PRB4] = ndgrid(DV_Ind,PR,PR,DV_Ind,PR,PR,DV_Ind,PR,PR,PR,PR); % Create a mesh grid of pressure to get all "combinations."

DV_Ind_HECM1 = DV_Ind_HECM1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA1 = PRA1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB1 = PRB1(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F1 = PRA1*ACap1-PRB1*ARod1;
MaxF1Diff = max(diff(sort(CPR_F1)))*RailFraction;

DV_Ind_HECM2 = DV_Ind_HECM2(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA2 = PRA2(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB2 = PRB2(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F2 = PRA2*ACap2-PRB2*ARod2;
MaxF2Diff = max(diff(sort(CPR_F2)))*RailFraction;

DV_Ind_HECM3 = DV_Ind_HECM3(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA3 = PRA3(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB3 = PRB3(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F3 = PRA3*ACap3-PRB3*ARod3;
MaxF3Diff = max(diff(sort(CPR_F3)))*RailFraction;

PRA4 = PRA4(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB4 = PRB4(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_T1 = (PRA4-PRB4)*D*ScaleHECM4;
MaxT1Diff = max(diff(sort(CPR_T1)))*RailFraction;

ones_t = ones(length(t),1);
ones_opt=ones(1,length(PRA1));
%% Find Angular Velocity of HECM Units @ Every Time
disp('HECM Hydraulics')

disp('HECM Hydraulics - P')

P1_Cyl = (ones_t*PRA1'*ACap1+F1*ones_opt)/ARod1;
DeltaP_HECM1 = P1_Cyl-ones_t*PRB1';
P2_Cyl = (ones_t*PRA2'*ACap2+F2*ones_opt)/ARod2;
DeltaP_HECM2 = P2_Cyl-ones_t*PRB2';
P3_Cyl = (ones_t*PRA3'*ACap3+F3*ones_opt)/ARod3;
DeltaP_HECM3 = P3_Cyl-ones_t*PRB3';
DeltaP_HECM4 = ones_t*(PRB4-PRA4)';  

disp('HECM Hydraulics - Q')

Q_HECM_1 = -repmat(V1,1,length(PRA1)).*ARod1;
Q_HECM_2 = -repmat(V2,1,length(PRA2)).*ARod2;
Q_HECM_3 = -repmat(V3,1,length(PRA3)).*ARod3;
Q_HECM_4 = ScaleHECM4*interp2(w_28cc_4Q, P_28cc_4Q, Q_28cc_4Q, W4*ones_opt,DeltaP_HECM4);

disp('HECM Hydraulics - w')

w1 = interp2(ScaleHECM1*Q_LinQ,P_LinQ,w_LinQ,Q_HECM_1,DeltaP_HECM1);
w2 = interp2(ScaleHECM2*Q_LinQ,P_LinQ,w_LinQ,Q_HECM_2,DeltaP_HECM2);
w3 = interp2(ScaleHECM3*Q_LinQ,P_LinQ,w_LinQ,Q_HECM_3,DeltaP_HECM3);
w4 = W4*ones_opt;

disp('HECM Hydraulics - T')
T1_Act = ScaleHECM1*interp2(w_28cc_4Q,P_28cc_4Q,T_28cc_4Q,w1,DeltaP_HECM1);
T2_Act = ScaleHECM2*interp2(w_28cc_4Q,P_28cc_4Q,T_28cc_4Q,w2,DeltaP_HECM2);
T3_Act = ScaleHECM3*interp2(w_28cc_4Q,P_28cc_4Q,T_28cc_4Q,w3,DeltaP_HECM3);
T4_Act = ScaleHECM4*interp2(w_28cc_4Q,P_28cc_4Q,T_28cc_4Q,W4*ones_opt,DeltaP_HECM4)-T4*ones_opt; 

%% Total HECM Hydraulic Power Loss
disp('HECM Hydraulics - Loss')

%     a1 = zeros(size(T1_Act));
%     a1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0) = DeltaP_HECM1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0).*Q_HECM_1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0);
%     a1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0) = -T1_Act(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0).*w1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0);
% 
%     b1 = T1_Act.*w1-DeltaP_HECM1.*Q_HECM_1;
%     b1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0) = T1_Act(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0).*w1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0);
%     b1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0) = -DeltaP_HECM1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0).*Q_HECM_1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0);
% 
%     PLoss_Hydraulic_HECM = b1-a1;
    
PLoss_Hydraulic_HECM = T1_Act.*w1-DeltaP_HECM1.*Q_HECM_1 + T2_Act.*w2-DeltaP_HECM2.*Q_HECM_2 + ...
    T3_Act.*w3-DeltaP_HECM3.*Q_HECM_3 + (T4*ones_opt+T4_Act).*w4-DeltaP_HECM4.*Q_HECM_4;

%% Electrical Losses HECM 1
disp('HECM electric')

load E_Map_ManikantaPalantla_3point2kW.mat

%comment out when there are no explicit torque limits set
PLoss_Electric_HECM1 = interp2(max(abs(w1(:)))/max(abs(w_m(:)))*w_m,MaxT1_Act/max(abs(T_m(:)))*T_m,max(abs(w1(:)))/max(abs(w_m(:)))*MaxT1_Act/max(abs(T_m(:)))*ELoss_m,w1,T1_Act);
PLoss_Electric_HECM2 = interp2(max(abs(w2(:)))/max(abs(w_m(:)))*w_m,MaxT2_Act/max(abs(T_m(:)))*T_m,max(abs(w2(:)))/max(abs(w_m(:)))*MaxT2_Act/max(abs(T_m(:)))*ELoss_m,w2,T2_Act);
PLoss_Electric_HECM3 = interp2(max(abs(w3(:)))/max(abs(w_m(:)))*w_m,MaxT3_Act/max(abs(T_m(:)))*T_m,max(abs(w3(:)))/max(abs(w_m(:)))*MaxT3_Act/max(abs(T_m(:)))*ELoss_m,w3,T3_Act);
PLoss_Electric_HECM4 = interp2(max(abs(w4(:)))/max(abs(w_m(:)))*w_m,MaxT4_Act/max(abs(T_m(:)))*T_m,max(abs(w4(:)))/max(abs(w_m(:)))*MaxT4_Act/max(abs(T_m(:)))*ELoss_m,w4,T4_Act);


%comment out when there are some torque limits mentioned:
% MaxT1_Act = max(abs(T1_Act(:))); %setting torque limits based on all possible torque values and on capped once
% MaxT2_Act = max(abs(T2_Act(:)));
% MaxT3_Act = max(abs(T3_Act(:)));
% MaxT4_Act = max(abs(T4_Act(:)));
% PLoss_Electric_HECM1 = interp2(max(abs(w1(:)))/max(abs(w_m(:)))*w_m,MaxT1_Act/max(abs(T_m(:)))*T_m,max(abs(w1(:)))/max(abs(w_m(:)))*MaxT1_Act/max(abs(T_m(:)))*ELoss_m,w1,T1_Act);
% PLoss_Electric_HECM2 = interp2(max(abs(w2(:)))/max(abs(w_m(:)))*w_m,MaxT2_Act/max(abs(T_m(:)))*T_m,max(abs(w2(:)))/max(abs(w_m(:)))*MaxT2_Act/max(abs(T_m(:)))*ELoss_m,w2,T2_Act);
% PLoss_Electric_HECM3 = interp2(max(abs(w3(:)))/max(abs(w_m(:)))*w_m,MaxT3_Act/max(abs(T_m(:)))*T_m,max(abs(w3(:)))/max(abs(w_m(:)))*MaxT3_Act/max(abs(T_m(:)))*ELoss_m,w3,T3_Act);
% PLoss_Electric_HECM4 = interp2(max(abs(w4(:)))/max(abs(w_m(:)))*w_m,MaxT4_Act/max(abs(T_m(:)))*T_m,max(abs(w4(:)))/max(abs(w_m(:)))*MaxT4_Act/max(abs(T_m(:)))*ELoss_m,w4,T4_Act);




PLoss_Electric_HECM = PLoss_Electric_HECM1 + PLoss_Electric_HECM2 + PLoss_Electric_HECM3 + PLoss_Electric_HECM4;

% battery_power = (T1_Act.*w1>0).*((T1_Act.*w1)+PLoss_Electric_HECM1) + (T1_Act.*w1<0).*((T1_Act.*w1)+PLoss_Electric_HECM1) + ...
%                 (T2_Act.*w2>0).*((T2_Act.*w2)+PLoss_Electric_HECM2) + (T2_Act.*w2<0).*((T2_Act.*w2)+PLoss_Electric_HECM2) + ...
%                 (T3_Act.*w3>0).*((T3_Act.*w3)+PLoss_Electric_HECM3) + (T3_Act.*w3<0).*((T3_Act.*w3)+PLoss_Electric_HECM3) + ...
%                 (T4_Act.*w4>0).*((T4_Act.*w4)+PLoss_Electric_HECM4) + (T4_Act.*w4<0).*((T4_Act.*w4)+PLoss_Electric_HECM4);
 
battery_power_HECM = (T1_Act.*w1>0).*((T1_Act.*w1)+PLoss_Electric_HECM1) + (T1_Act.*w1<0).*((T1_Act.*w1)+PLoss_Electric_HECM1) + ...
                (T2_Act.*w2>0).*((T2_Act.*w2)+PLoss_Electric_HECM2) + (T2_Act.*w2<0).*((T2_Act.*w2)+PLoss_Electric_HECM2) + ...
                (T3_Act.*w3>0).*((T3_Act.*w3)+PLoss_Electric_HECM3) + (T3_Act.*w3<0).*((T3_Act.*w3)+PLoss_Electric_HECM3) + ...
                (T4_Act.*w4>0).*((T4_Act.*w4)+PLoss_Electric_HECM4) + (T4_Act.*w4<0).*((T4_Act.*w4)+PLoss_Electric_HECM4);
 

% ElectricEff1 = 0.9; ElectricEff2 = 0.9;
% ElectricEff3 = 0.9; ElectricEff4 = 0.9;
% 
% battery_power = (T1_Act.*w1>0).*(T1_Act.*w1)/ElectricEff1+(T1_Act.*w1<0).*(T1_Act.*w1)*ElectricEff1 + ...
%                 (T2_Act.*w2>0).*(T2_Act.*w2)/ElectricEff2+(T2_Act.*w2<0).*(T2_Act.*w2)*ElectricEff2 + ...
%                 (T3_Act.*w3>0).*(T3_Act.*w3)/ElectricEff3+(T3_Act.*w3<0).*(T3_Act.*w3)*ElectricEff3 + ...
%                 (T4_Act.*w4>0).*(T4_Act.*w4)/ElectricEff4+(T4_Act.*w4<0).*(T4_Act.*w4)*ElectricEff4;
% 
% PLoss_Electric_HECM = (T1_Act.*w1>0).*abs(T1_Act.*w1*(1/ElectricEff1-1))...
%                      +(T1_Act.*w1<0).*abs(T1_Act.*w1*(1-ElectricEff1))+...
%                      (T2_Act.*w2>0).*abs(T2_Act.*w2*(1/ElectricEff2-1))...
%                      +(T2_Act.*w2<0).*abs(T2_Act.*w2*(1-ElectricEff2))+...
%                      (T3_Act.*w3>0).*abs(T3_Act.*w3*(1/ElectricEff3-1))...
%                      +(T3_Act.*w3<0).*abs(T3_Act.*w3*(1-ElectricEff3))+...
%                      (T4_Act.*w4>0).*abs(T4_Act.*w4*(1/ElectricEff4-1))...
%                      +(T4_Act.*w4<0).*abs(T4_Act.*w4*(1-ElectricEff4));

%% Find Main Pump Flows and Losses (Rail Losses)
disp('Main pump losses');

load AP_107cc_Sept16.mat

wMP = 2000*2*pi/60; % Angular Velocity [rad/s]
T_PR = interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
% T_PR = interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,[PR';-PR'],ones(2*(length(PR)),1)*wMP);
Q_PR = interp2(P1_Mapping,w1rad_Mapping,Q1_Act_Mapping,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
% Q_PR = interp2(P1_Mapping,w1rad_Mapping,Q1_Act_Mapping,[PR';-PR'],ones(2*(length(PR)),1)*wMP);
%% Finding R values involving both the losses:
MP_Hyd_Loss = T_PR*wMP - Q_PR.*[PR(2:end)'; -PR(2:end)'];
% MP_Hyd_Loss = T_PR*wMP - Q_PR.*[PR'; -PR'];
MP_Elec_Loss = interp2((wMP/max(abs(w_m(:))))*w_m, (max(abs(T_PR))/max(abs(T_m(:))))*T_m, ((wMP/max(abs(w_m(:))))*(max(abs(T_PR))/max(abs(T_m(:))))*ELoss_m), wMP*ones(size(T_PR)), T_PR);

%calculating the efficiency of the main pump/generator for the pumping and the motoring case: 
Pump_Eff = (Q_PR(1:length(PR)-1).*[PR(2:end)'])./((Q_PR(1:length(PR)-1).*[PR(2:end)']) + MP_Hyd_Loss(1:length(PR)-1) + MP_Elec_Loss(1:length(PR)-1));
% Pump_Eff = (Q_PR(1:length(PR)).*[PR'])./((Q_PR(1:length(PR)).*[PR']) + MP_Hyd_Loss(1:length(PR)) + MP_Elec_Loss(1:length(PR)));
Motor_Eff = (-Q_PR(length(PR):end).*[-PR(2:end)'] - MP_Hyd_Loss(length(PR):end) - MP_Elec_Loss(length(PR):end))./(-Q_PR(length(PR):end).*[-PR(2:end)']);
% Motor_Eff = (-Q_PR((length(PR)+1):end).*[-PR'] - MP_Hyd_Loss((length(PR)+1):end) - MP_Elec_Loss((length(PR)+1):end))./(-Q_PR((length(PR)+1):end).*[-PR']);

%Using the efficiency values to calculate R values. Per unit power:
R_P_loss = (1./Pump_Eff - 1);
R_M_loss = (1 - Motor_Eff);

%  MP_PumpEff = (Q_PR(1:(length(PR)-1)).*PR(2:end)')./(T_PR(1:(length(PR)-1))*wMP);
%  MP_MotorEff = T_PR(length(PR):end)*wMP./(-PR(2:end)'.*Q_PR(length(PR):end));

%MP_loss=(T_PR*wMP-Q_PR.*[PR(2:3)';-PR(2:3)'])./(Q_PR.*[PR(2:3)';PR(2:3)'])

% R_P_loss = (1./MP_PumpEff-1);  % loss per power
% R_M_loss = (1-MP_MotorEff);

% R1_P_loss = (1/MP_PumpEff(1)-1);  % loss per power
% R1_M_loss = (1-MP_MotorEff(1));
% R2_P_loss = (1/MP_PumpEff(2)-1);
% R2_M_loss = (1-MP_MotorEff(2));
% if length(PR)>3,
%     R3_P_loss = (1/MP_PumpEff(3)-1);
%     R3_M_loss = (1-MP_MotorEff(3));
% end

%QR1 = zeros(length(t),length(PRA1));
%QR2 = zeros(length(t),length(PRA1));
%if length(PR)>3,
%    QR3 = zeros(length(t),length(PRA1));
%end

%finding flows for each of the rails
for k=1:length(PR),
    QR{k} = zeros(length(t),length(PRA1));
    ind = find(PRA1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+V1*ACap1*ones(1,length(ind));
    ind = find(PRB1==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+Q_HECM_1(:,ind);
    ind = find(PRA2==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+V2*ACap2*ones(1,length(ind));
    ind = find(PRB2==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+Q_HECM_2(:,ind);
    ind = find(PRA3==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+V3*ACap3*ones(1,length(ind));
    ind = find(PRB3==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+Q_HECM_3(:,ind);
    ind = find(PRA4==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)+Q_HECM_4(:,ind);
    ind = find(PRB4==PR(k)); QR{k}(:,ind) = QR{k}(:,ind)-Q_HECM_4(:,ind);
end

% ind = find(PRA1==PR(2)); QR1(:,ind) = QR1(:,ind)+V1*ACap1*ones(1,length(ind));
% ind = find(PRB1==PR(2)); QR1(:,ind) = QR1(:,ind)+Q_HECM_1(:,ind);
% ind = find(PRA2==PR(2)); QR1(:,ind) = QR1(:,ind)+V2*ACap2*ones(1,length(ind));
% ind = find(PRB2==PR(2)); QR1(:,ind) = QR1(:,ind)+Q_HECM_2(:,ind);
% 
% ind = find(PRA1==PR(3)); QR2(:,ind) = QR2(:,ind)+V1*ACap1*ones(1,length(ind));
% ind = find(PRB1==PR(3)); QR2(:,ind) = QR2(:,ind)+Q_HECM_1(:,ind);
% ind = find(PRA2==PR(3)); QR2(:,ind) = QR2(:,ind)+V2*ACap2*ones(1,length(ind));
% ind = find(PRB2==PR(3)); QR2(:,ind) = QR2(:,ind)+Q_HECM_2(:,ind);
% 
% ind = find(PRA1==PR(4)); QR3(:,ind) = QR3(:,ind)+V1*ACap1*ones(1,length(ind));
% ind = find(PRB1==PR(4)); QR3(:,ind) = QR3(:,ind)+Q_HECM_1(:,ind);
% ind = find(PRA2==PR(4)); QR3(:,ind) = QR3(:,ind)+V2*ACap2*ones(1,length(ind));
% ind = find(PRB2==PR(4)); QR3(:,ind) = QR3(:,ind)+Q_HECM_2(:,ind);


%% battery power corresponding to the main pump:

battery_power_rail_p = cell(1,length(PR)-1);
% battery_power_rail_p = cell(1,length(PR));
battery_power_rail_m = cell(1,length(PR)-1);
% battery_power_rail_m = cell(1,length(PR));

for i = 2:length(PR)
    battery_power_rail_p{i-1} = QR{i}*PR(i)*(1+R_P_loss(i-1));
    battery_power_rail_m{i-1} = QR{i}*PR(i)*(1-R_M_loss(i-1));
end

% for i = 1:length(PR)
%     battery_power_rail_p{i} = QR{i}*PR(i)*(1+R_P_loss(i));
%     battery_power_rail_m{i} = QR{i}*PR(i)*(1-R_M_loss(i));
% end

% battery_power_rail_net = zeros(size(battery_power_HECM));
% for i = 1:length(battery_power_rail)
%     battery_power_rail_net = battery_power_rail_net + battery_power_rail{i};
% end

% %% net battery power
% battery_power = battery_power_HECM + battery_power_rail_net;

%% Sum All Losses/Costs ***************************************************
disp('Constraints')
% Torque limit

% %comment out if max limits are set already
PLoss_Electric_HECM(abs(T1_Act) > MaxT1_Act | abs(T2_Act) > MaxT2_Act | abs(T3_Act) > MaxT3_Act | abs(T4_Act) > MaxT4_Act) = inf;


%PLoss_Electric_HECM(abs(T4_Act) > MaxT4_Act) = inf;
% Cavitation limit
PLoss_Hydraulic_HECM(P1_Cyl < -1e5 | P2_Cyl < -1e5 | P3_Cyl < -1e5) = inf;
HECMLosses = PLoss_Hydraulic_HECM + PLoss_Electric_HECM;
% battery_power; 
%QR1*PR(2)*R1_P_loss + QR2*PR(3)*R2_M_loss
%%
% toc
