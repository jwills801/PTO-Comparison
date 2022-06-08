% function [TotalCost, I] = InnerLoopOptFunction_MultHECMs(lambda1)
% Jake's original HHEA script that I adapated for wec-sim
tic
load JCB_LongerDC.mat
% load jcb48ZI.mat

startingInd = 1;
spacing = 8;
ScaleHECM1 = .14;
ScaleHECM2 = .21;
ScaleHECM3 = .18;
ScaleHECM4 = 18.1/107;

% lambda1 =    0.249256134033202; % One Pressure Rail + RodSideONLY + 90 Deg
% lambda1 =    -0.063618397712708; % One Pressure Rail + RodSideONLY + Grading
% lambda1 =    0.271308708190917; % One Pressure Rail + RodSideONLY + Trenching

% lambda1 =       0.243463802337646; % One Pressure Rail + 90 Deg
% lambda1 =       0.073869228363036; % One Pressure Rail + Grading
lambda1 =         0.244643402099609; % One Pressure Rail + Trenching

maxTorque = 350;

SwingGearRatio = Swing_ratio;

%% Test 4 HECM Unit

%% Fluid Properties
mu = (32e-6)*870; % Dynamic Viscosity []
B = 1.7e9; % Bulk Modulus []
rho = 870; % Density [kg/m^3]

%% Main Pump Properties

[Pmap_MP, Qmap_MP, Wmap_MP, DMP, CfMP, ChMP, CvMP, CsMP, CstMP] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');

xMP=1; % Frac Disp
wMP = 2000*2*pi/60; % Angular Velocity [rad/s]
%% HECM Pump 1 Properties

[Pmap1, Qmap1, Wmap1, D1, Cf1, Ch1, Cv1, Cs1, Cst1,] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');
x1=1; % Frac Disp

%% HECM Pump 2 Properties

[Pmap2, Qmap2, Wmap2, D2, Cf2, Ch2, Cv2, Cs2, Cst2] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');
x2=1; % Frac Disp

%% HECM Pump 3 Properties

[Pmap3, Qmap3, Wmap3, D3, Cf3, Ch3, Cv3, Cs3, Cst3] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');
x3=1; % Frac Disp

%% HECM Pump 4 Properties

[Pmap4, Qmap4, Wmap4, D4, Cf4, Ch4, Cv4, Cs4, Cst4] = PumpInfo('AP_107cc_M_FullDisp_TESTpumpINFOfunction');
x4=1; % Frac Disp

%%
load('AP_107cc_M_FullDisp_TESTpumpINFOfunction.mat')

%% Load Drive Cycle/ Actuator Dimensions / HECM Dimensions

% load JCB Small Excavator Drive Cycle
% *** NOTE that the forces are flipped. This is because the data file
% provided forces with a positive sign when acting in the direction from
% rod to cap end. The following code assumed a positive force when acting
% from cap to rod end on the cylinder.
load JCB_LongerDC.mat
% load jcb48ZI.mat

% Key:
% LinearActuator 1: Boom
% LinearActuator 2: Arm
% LinearActuator 3: Bucket
% RotaryActuator 1: Swing

% % 90 Degree Drive Cycle***************************************************
% t_startInd = 1;
% t_endInd = length(NinetyDeg.time);
% 
% t = NinetyDeg.time(t_startInd:spacing:t_endInd);
% t_step = t(2)-t(1);
% 
% F1 = -NinetyDeg.boom_f(t_startInd:spacing:t_endInd);
% V1 = NinetyDeg.boom_v(t_startInd:spacing:t_endInd);
% ACap1 = boom_Ac;
% ARod1 = boom_Ar;
% 
% F2 = -NinetyDeg.arm_f(t_startInd:spacing:t_endInd);
% V2 = NinetyDeg.arm_v(t_startInd:spacing:t_endInd);
% ACap2 = arm_Ac;
% ARod2 = arm_Ar;
% 
% F3 = -NinetyDeg.bkt_f(t_startInd:spacing:t_endInd);
% V3 = NinetyDeg.bkt_v(t_startInd:spacing:t_endInd);
% ACap3 = bkt_Ac;
% ARod3 = bkt_Ar;
% 
% T1 = -NinetyDeg.sw_t(t_startInd:spacing:t_endInd)/SwingGearRatio;
% W1 = NinetyDeg.sw_w(t_startInd:spacing:t_endInd)*SwingGearRatio;



% % Grading Drive Cycle *****************************************************
% t_startInd = 1;
% t_endInd = length(Grading.time);
% 
% t = Grading.time(t_startInd:spacing:t_endInd);
% t_step = t(2)-t(1);
% 
% F1 = -Grading.boom_f(t_startInd:spacing:t_endInd);
% V1 = Grading.boom_v(t_startInd:spacing:t_endInd);
% ACap1 = boom_Ac;
% ARod1 = boom_Ar;
% 
% F2 = -Grading.arm_f(t_startInd:spacing:t_endInd);
% V2 = Grading.arm_v(t_startInd:spacing:t_endInd);
% ACap2 = arm_Ac;
% ARod2 = arm_Ar;
% 
% F3 = -Grading.bkt_f(t_startInd:spacing:t_endInd);
% V3 = Grading.bkt_v(t_startInd:spacing:t_endInd);
% ACap3 = bkt_Ac;
% ARod3 = bkt_Ar;
% 
% T1 = -Grading.sw_t(t_startInd:spacing:t_endInd)/SwingGearRatio;
% W1 = Grading.sw_w(t_startInd:spacing:t_endInd)*SwingGearRatio;



% Trenching Drive Cycle *****************************************************
t_startInd = 1;
t_endInd = length(Trenching.time);

t = Trenching.time(t_startInd:spacing:t_endInd);
t_step = t(2)-t(1);

F1 = -Trenching.boom_f(t_startInd:spacing:t_endInd);
V1 = Trenching.boom_v(t_startInd:spacing:t_endInd);
ACap1 = boom_Ac;
ARod1 = boom_Ar;

F2 = -Trenching.arm_f(t_startInd:spacing:t_endInd);
V2 = Trenching.arm_v(t_startInd:spacing:t_endInd);
ACap2 = arm_Ac;
ARod2 = arm_Ar;

F3 = -Trenching.bkt_f(t_startInd:spacing:t_endInd);
V3 = Trenching.bkt_v(t_startInd:spacing:t_endInd);
ACap3 = bkt_Ac;
ARod3 = bkt_Ar;

T1 = -Trenching.sw_t(t_startInd:spacing:t_endInd)/SwingGearRatio;
W1 = Trenching.sw_w(t_startInd:spacing:t_endInd)*SwingGearRatio;

%% Given PR Selection and Directional Valve Indicator

% PR = [0 1/2 1]*29e6; % Define an array of pressures
PR = [0 1]*29e6; % Define an array of pressures
DV_Ind = [0 1];

%% Find Pressure Differentials For HECMs, Pressures in the Cylinders (Depending on Which Combination of Rails/Valve and the Force)

[DV_Ind_HECM1, PRA1, PRB1, DV_Ind_HECM2, PRA2, PRB2, DV_Ind_HECM3, PRA3, PRB3, PRA4, PRB4] = ndgrid(DV_Ind,PR,PR,DV_Ind,PR,PR,DV_Ind,PR,PR,PR,PR); % Create a mesh grid of pressure to get all "combinations."

% IndexMatrix = zeros(size(DV_Ind_HECM1));

DV_Ind_HECM1 = DV_Ind_HECM1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA1 = PRA1(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB1 = PRB1(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F1 = PRA1*ACap1-PRB1*ARod1;

DV_Ind_HECM2 = DV_Ind_HECM2(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA2 = PRA2(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB2 = PRB2(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F2 = PRA2*ACap2-PRB2*ARod2;

DV_Ind_HECM3 = DV_Ind_HECM3(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRA3 = PRA3(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB3 = PRB3(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_F3 = PRA3*ACap3-PRB3*ARod3;

PRA4 = PRA4(:); % Spread the combinations of pressure rails and directional valves into a vector.
PRB4 = PRB4(:); % Spread the combinations of pressure rails and directional valves into a vector.
CPR_T1 = (PRA4-PRB4)*D4*ScaleHECM4;

PA1 = PRA1; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on B Side)
PB1 = PRB1; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on A Side)

PA2 = PRA2; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on B Side)
PB2 = PRB2; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on A Side)

PA3 = PRA3; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on B Side)
PB3 = PRB3; % Creates an vector of pressures for the side of the actuator not connected to HECM. (Assumes HECM is on A Side)

CPR_A1 = repmat(PA1',length(t),1);
CPR_B1 = repmat(PB1',length(t),1);
CPR_A2 = repmat(PA2',length(t),1);
CPR_B2 = repmat(PB2',length(t),1);
CPR_A3 = repmat(PA3',length(t),1);
CPR_B3 = repmat(PB3',length(t),1);
CPR_A4 = repmat(PRA4',length(t),1);
CPR_B4 = repmat(PRB4',length(t),1);


P1_Cyl = zeros(length(t),length(PRA1));
P2_Cyl = zeros(length(t),length(PRA1));
P3_Cyl = zeros(length(t),length(PRA1));
DeltaP_HECM1 = zeros(length(t),length(PRA1));
DeltaP_HECM2 = zeros(length(t),length(PRA1));
DeltaP_HECM3 = zeros(length(t),length(PRA1));

for i = 1:size(PRA1,1)
    
    if DV_Ind_HECM1(i) == 0
        P1_Cyl(:,i) = (PB1(i)*ARod1-F1)/ACap1;
        DeltaP_HECM1(:,i) = P1_Cyl(:,i)-PRA1(i);
    else
        P1_Cyl(:,i) = (PA1(i)*ACap1+F1)/ARod1;
        DeltaP_HECM1(:,i) = P1_Cyl(:,i)-PRB1(i);
    end
    
    
    if DV_Ind_HECM2(i) == 0
        P2_Cyl(:,i) = (PB2(i)*ARod2-F2)/ACap2;
        DeltaP_HECM2(:,i) = P2_Cyl(:,i)-PRA2(i);
    else
        P2_Cyl(:,i) = (PA2(i)*ACap2+F2)/ARod2;
        DeltaP_HECM2(:,i) = P2_Cyl(:,i)-PRB2(i);
    end
        
    if DV_Ind_HECM3(i) == 0
        P3_Cyl(:,i) = (PB3(i)*ARod3-F3)/ACap3;
        DeltaP_HECM3(:,i) = P3_Cyl(:,i)-PRA3(i);
    else
        P3_Cyl(:,i) = (PA3(i)*ACap3+F3)/ARod3;
        DeltaP_HECM3(:,i) = P3_Cyl(:,i)-PRB3(i);
    end
    
end

% DeltaP_HECM4 = CPR_A4-CPR_B4;
DeltaP_HECM4 = CPR_B4-CPR_A4;

%% Find Angular Velocity of HECM Units @ Every Time
V1_Array = repmat(V1,1,length(PRA1)); % Repeat the desired velocity matrix for next calculation
Q1_A = V1_Array.*ACap1; % The flow taken from rail connected to the A side.(HECM1)
Q1_B = -V1_Array.*ARod1; % The flow taken from rail connected to the B side.(HECM1)

V2_Array = repmat(V2,1,length(PRA1)); % Repeat the desired velocity matrix for next calculation
Q2_A = V2_Array.*ACap2; % The flow taken from rail connected to the A side.(HECM2)
Q2_B = -V2_Array.*ARod2; % The flow taken from rail connected to the B side.(HECM2)

V3_Array = repmat(V3,1,length(PRA1)); % Repeat the desired velocity matrix for next calculation
Q3_A = V3_Array.*ACap3; % The flow taken from rail connected to the A side.(HECM3)
Q3_B = -V3_Array.*ARod3; % The flow taken from rail connected to the B side.(HECM4)

Q_HECM_1 = zeros(length(t),length(PRA1));
Q_HECM_2 = zeros(length(t),length(PRA1));
Q_HECM_3 = zeros(length(t),length(PRA1));

for i = 1:length(DV_Ind_HECM1)
    
    if DV_Ind_HECM1(i) == 0
        Q_HECM_1(:,i) = Q1_A(:,i);
    else
        Q_HECM_1(:,i) = Q1_B(:,i);
    end
    
    if DV_Ind_HECM2(i) == 0
        Q_HECM_2(:,i) = Q2_A(:,i);
    else
        Q_HECM_2(:,i) = Q2_B(:,i);
    end

    if DV_Ind_HECM3(i) == 0
        Q_HECM_3(:,i) = Q3_A(:,i);
    else
        Q_HECM_3(:,i) = Q3_B(:,i);
    end

end

w1 = interp2(Pmap1,ScaleHECM1*Qmap1,Wmap1,DeltaP_HECM1,Q_HECM_1);
w2 = interp2(Pmap2,ScaleHECM2*Qmap2,Wmap2,DeltaP_HECM2,Q_HECM_2);
w3 = interp2(Pmap3,ScaleHECM3*Qmap3,Wmap3,DeltaP_HECM3,Q_HECM_3);
w4 = repmat(W1,1,size(w3,2));

Q_HECM_4 = ScaleHECM4*interp2(P1_Mapping,w1rad_Mapping,Q1_Act_Mapping,DeltaP_HECM4,w4);
    
T1_Ideal = DeltaP_HECM1.*D1*ScaleHECM1;
T2_Ideal = DeltaP_HECM2.*D2*ScaleHECM2;
T3_Ideal = DeltaP_HECM3.*D3*ScaleHECM3;
T4_Ideal = DeltaP_HECM4.*D4*ScaleHECM4;

T1_Act = ScaleHECM1*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM1,w1);
T2_Act = ScaleHECM2*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM2,w2);
T3_Act = ScaleHECM3*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM3,w3);
T4_Act = ScaleHECM4*interp2(P1_Mapping,w1rad_Mapping,T1_Act_Mapping,DeltaP_HECM4,w4)-repmat(T1,1,size(T4_Ideal,2));

%% Find HECM 1 Pump Losses

% Find the flow and torque loss for the HECM operating at nominal
% conditions.
% Loss1_Q = ScaleHECM1*(((abs(D1*Cs1*(DeltaP_HECM1)/mu) + abs(D1.*w1.*x1.*(DeltaP_HECM1)./B) + D1^(2/3)*Cst1*(2*abs(DeltaP_HECM1)/rho).^.5)));
% Loss1_T = ScaleHECM1*(((abs(D1*Cv1*mu*w1) + abs(D1*(DeltaP_HECM1)*Cf1) + abs(Ch1*x1^3*w1.^2*rho*D1^(5/3)/2))));
Loss1_Q = abs(Q_HECM_1-w1*D1*ScaleHECM1);
Loss1_T = abs(T1_Act-T1_Ideal);

% Find the power loss of the HECM.
PLoss1_Q = abs(DeltaP_HECM1).*Loss1_Q;
PLoss1_T = abs(w1).*Loss1_T;
PLoss1_Hydraulic = PLoss1_Q+PLoss1_T;

%% Find HECM 2 Pump Losses

% Find the flow and torque loss for the HECM operating at nominal
% conditions.
% Loss2_Q = ScaleHECM2*(((abs(D2*Cs2*(DeltaP_HECM2)/mu) + abs(D2.*w2.*x2.*(DeltaP_HECM2)./B) + D1^(2/3)*Cst2*(2*abs(DeltaP_HECM2)/rho).^.5)));
% Loss2_T = ScaleHECM2*(((abs(D2*Cv2*mu*w2) + abs(D2*(DeltaP_HECM2)*Cf2) + abs(Ch2*x2^3*w2.^2*rho*D2^(5/3)/2))));
Loss2_Q = abs(Q_HECM_2-w2*D2*ScaleHECM2);
Loss2_T = abs(T2_Act-T2_Ideal);

% Find the power loss of the HECM.
PLoss2_Q = abs(DeltaP_HECM2).*Loss2_Q;
PLoss2_T = abs(w2).*Loss2_T;
PLoss2_Hydraulic = PLoss2_Q+PLoss2_T;

%% Find HECM 3 Pump Losses

% Find the flow and torque loss for the HECM operating at nominal
% conditions.
% Loss3_Q = ScaleHECM3*(((abs(D3*Cs3*(DeltaP_HECM3)/mu) + abs(D3.*w3.*x3.*(DeltaP_HECM3)./B) + D3^(2/3)*Cst3*(2*abs(DeltaP_HECM3)/rho).^.5)));
% Loss3_T = ScaleHECM3*(((abs(D3*Cv3*mu*w3) + abs(D3*(DeltaP_HECM3)*Cf3) + abs(Ch3*x3^3*w3.^2*rho*D3^(5/3)/2))));
Loss3_Q = abs(Q_HECM_3-w3*D3*ScaleHECM3);
Loss3_T = abs(T3_Act-T3_Ideal);

% Find the power loss of the HECM.
PLoss3_Q = abs(DeltaP_HECM3).*Loss3_Q;
PLoss3_T = abs(w3).*Loss3_T;
PLoss3_Hydraulic = PLoss3_Q+PLoss3_T;

%% Find HECM 4 Pump Losses

% Find the flow and torque loss for the HECM operating at nominal
% conditions.
% Loss4_Q = ScaleHECM4*(((abs(D4*Cs4*(DeltaP_HECM4)/mu) + abs(D4.*w4.*x4.*(DeltaP_HECM4)./B) + D4^(2/3)*Cst4*(2*abs(DeltaP_HECM4)/rho).^.5)));
% Loss4_T = ScaleHECM4*(((abs(D4*Cv4*mu*w4) + abs(D4*(DeltaP_HECM4)*Cf4) + abs(Ch4*x4^3*w4.^2*rho*D4^(5/3)/2))));
Loss4_Q = abs(Q_HECM_4-w4*D4*ScaleHECM4);
Loss4_T = abs(T4_Act-T4_Ideal-T1);

% Find the power loss of the HECM.
PLoss4_Q = abs(DeltaP_HECM4).*Loss4_Q;
PLoss4_T = abs(w4).*Loss4_T;
PLoss4_Hydraulic = PLoss4_Q+PLoss4_T;

% Find the Actual Flow Rate Across the Rotary HECM
Q4_A = Q_HECM_4;
Q4_B = -Q4_A;

%% Total HECM Hydraulic Power Loss

% PLoss_Hydraulic_HECM = PLoss1_Hydraulic+PLoss2_Hydraulic+PLoss3_Hydraulic+PLoss4_Hydraulic;

PLoss_Hydraulic_HECM1_check = abs(T1_Act.*w1-DeltaP_HECM1.*Q_HECM_1);
PLoss_Hydraulic_HECM2_check = abs(T2_Act.*w2-DeltaP_HECM2.*Q_HECM_2);
PLoss_Hydraulic_HECM3_check = abs(T3_Act.*w3-DeltaP_HECM3.*Q_HECM_3);
PLoss_Hydraulic_HECM4_check = abs((T4_Act+repmat(T1,1,size(T4_Ideal,2))).*w4-DeltaP_HECM4.*Q_HECM_4);

PLoss_Hydraulic_HECM = PLoss_Hydraulic_HECM1_check+PLoss_Hydraulic_HECM2_check+PLoss_Hydraulic_HECM3_check+PLoss_Hydraulic_HECM4_check;

%% Electrical Losses HECM 1
% Assume a constant electrical efficiency.
ElectricEff1 = .9;

% The following process will allow me to sum the power losses when pumping
% and motoring and not "double count" anything. It is necessary to consider
% these cases (motoring/pumping) independently as the losses are calculated
% differently for either case.

w1_Generator = w1; % Grabs the angular velocity of HECM
w1_Generator(sign(T1_Act.*w1)>=0) = 0; % Throws away values when Motoring
w1_Motor = w1; % Grabs the angular velocity of HECM
w1_Motor(sign(T1_Act.*w1)<0) = 0; % Throws away values when Generating

PLoss1_Electric_Generator = abs(T1_Act.*w1_Generator*(1-ElectricEff1));
PLoss1_Electric_Motor = abs(T1_Act.*w1_Motor*(1/ElectricEff1-1));

PLoss1_Electric = PLoss1_Electric_Generator+PLoss1_Electric_Motor;

%% Electrical Losses HECM 2
% Assume a constant electrical efficiency.
ElectricEff2 = .9;

% The following process will allow me to sum the power losses when pumping
% and motoring and not "double count" anything. It is necessary to consider
% these cases (motoring/pumping) independently as the losses are calculated
% differently for either case.

w2_Generator = w2; % Grabs the angular velocity of HECM
w2_Generator(sign(T2_Act.*w2)>=0) = 0; % Throws away values when Motoring
w2_Motor = w2; % Grabs the angular velocity of HECM
w2_Motor(sign(T2_Act.*w2)<0) = 0; % Throws away values when Generating

PLoss2_Electric_Generator = abs(T2_Act.*w2_Generator*(1-ElectricEff2));
PLoss2_Electric_Motor = abs(T2_Act.*w2_Motor*(1/ElectricEff2-1));

PLoss2_Electric = PLoss2_Electric_Generator+PLoss2_Electric_Motor;

%% Electrical Losses HECM 3
% Assume a constant electrical efficiency.
ElectricEff3 = .9;

% The following process will allow me to sum the power losses when pumping
% and motoring and not "double count" anything. It is necessary to consider
% these cases (motoring/pumping) independently as the losses are calculated
% differently for either case.

w3_Generator = w3; % Grabs the angular velocity of HECM
w3_Generator(sign(T3_Act.*w3)>=0) = 0; % Throws away values when Motoring
w3_Motor = w3; % Grabs the angular velocity of HECM
w3_Motor(sign(T3_Act.*w3)<0) = 0; % Throws away values when Generating

PLoss3_Electric_Generator = abs(T3_Act.*w3_Generator*(1-ElectricEff3));
PLoss3_Electric_Motor = abs(T3_Act.*w3_Motor*(1/ElectricEff3-1));

PLoss3_Electric = PLoss3_Electric_Generator+PLoss3_Electric_Motor;

%% Electrical Losses HECM 4
% Assume a constant electrical efficiency.
ElectricEff4 = .9;

% The following process will allow me to sum the power losses when pumping
% and motoring and not "double count" anything. It is necessary to consider
% these cases (motoring/pumping) independently as the losses are calculated
% differently for either case.

w4_Generator = w4; % Grabs the angular velocity of HECM
w4_Generator(sign(T4_Act.*w4)>=0) = 0; % Throws away values when Motoring
w4_Motor = w4; % Grabs the angular velocity of HECM
w4_Motor(sign(T4_Act.*w4)<0) = 0; % Throws away values when Generating

PLoss4_Electric_Generator = abs(T4_Act.*w4_Generator*(1-ElectricEff4));
PLoss4_Electric_Motor = abs(T4_Act.*w4_Motor*(1/ElectricEff4-1));

PLoss4_Electric = PLoss4_Electric_Generator+PLoss4_Electric_Motor;
PLoss4_Electric(abs(T4_Act)>maxTorque)=inf;

%% Total Electric Power Loss

PLoss_Electric_HECM = PLoss1_Electric+PLoss2_Electric+PLoss3_Electric+PLoss4_Electric;

% %% Find Main Pump Losses (Rail Losses)
% 
% % Find the flow and torque loss for the main pump operating at nominal
% % conditions. HECM1
% LossMP_QA1 = ((abs(DMP*CsMP*(CPR_A1)/mu) + abs(DMP.*wMP.*xMP.*(CPR_A1)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_A1)/rho).^.5)); 
% LossMP_TA1 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_A1)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% LossMP_QB1 = ((abs(DMP*CsMP*(CPR_B1)/mu) + abs(DMP.*wMP.*xMP.*(CPR_B1)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_B1)/rho).^.5));
% LossMP_TB1 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_B1)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% 
% % Find the flow and torque loss for the main pump operating at nominal
% % conditions. HECM2
% LossMP_QA2 = ((abs(DMP*CsMP*(CPR_A2)/mu) + abs(DMP.*wMP.*xMP.*(CPR_A2)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_A2)/rho).^.5)); 
% LossMP_TA2 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_A2)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% LossMP_QB2 = ((abs(DMP*CsMP*(CPR_B2)/mu) + abs(DMP.*wMP.*xMP.*(CPR_B2)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_B2)/rho).^.5));
% LossMP_TB2 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_B2)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% 
% % Find the flow and torque loss for the main pump operating at nominal
% % conditions. HECM3
% LossMP_QA3 = ((abs(DMP*CsMP*(CPR_A3)/mu) + abs(DMP.*wMP.*xMP.*(CPR_A3)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_A3)/rho).^.5)); 
% LossMP_TA3 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_A3)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% LossMP_QB3 = ((abs(DMP*CsMP*(CPR_B3)/mu) + abs(DMP.*wMP.*xMP.*(CPR_B3)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_B3)/rho).^.5));
% LossMP_TB3 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_B3)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% 
% % Find the flow and torque loss for the main pump operating at nominal
% % conditions. HECM4 (Torque 1)
% LossMP_QA4 = ((abs(DMP*CsMP*(CPR_A4)/mu) + abs(DMP.*wMP.*xMP.*(CPR_A4)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_A4)/rho).^.5)); 
% LossMP_TA4 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_A4)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% LossMP_QB4 = ((abs(DMP*CsMP*(CPR_B4)/mu) + abs(DMP.*wMP.*xMP.*(CPR_B4)./B) + DMP^(2/3)*CstMP*(2*abs(CPR_B4)/rho).^.5));
% LossMP_TB4 = ((abs(DMP*CvMP*mu*wMP) + abs(DMP*(CPR_B4)*CfMP) + abs(ChMP*xMP^3*wMP.^2*rho*DMP^(5/3)/2)));
% 
% % Find the power loss of the main pump. HECM1 (& set Tank Losses to 0)
% PLossMP_A1 = LossMP_QA1.*abs(CPR_A1)+LossMP_TA1*abs(wMP);
% PLossMP_B1 = LossMP_QB1.*abs(CPR_B1)+LossMP_TB1*abs(wMP);
% PLossMP_A1(CPR_A1==0)=0;
% PLossMP_B1(CPR_B1==0)=0;
% 
% % Find the power loss of the main pump. HECM2 (& set Tank Losses to 0)
% PLossMP_A2 = LossMP_QA2.*abs(CPR_A2)+LossMP_TA2*abs(wMP);
% PLossMP_B2 = LossMP_QB2.*abs(CPR_B2)+LossMP_TB2*abs(wMP);
% PLossMP_A2(CPR_A2==0)=0;
% PLossMP_B2(CPR_B2==0)=0;
% 
% % Find the power loss of the main pump. HECM3 (& set Tank Losses to 0)
% PLossMP_A3 = LossMP_QA3.*abs(CPR_A3)+LossMP_TA3*abs(wMP);
% PLossMP_B3 = LossMP_QB3.*abs(CPR_B3)+LossMP_TB3*abs(wMP);
% PLossMP_A3(CPR_A3==0)=0;
% PLossMP_B3(CPR_B3==0)=0;
% 
% % Find the power loss of the main pump. HECM4 (& set Tank Losses to 0)
% PLossMP_A4 = LossMP_QA4.*abs(CPR_A4)+LossMP_TA4*abs(wMP);
% PLossMP_B4 = LossMP_QB4.*abs(CPR_B4)+LossMP_TB4*abs(wMP);
% PLossMP_A4(CPR_A4==0)=0;
% PLossMP_B4(CPR_B4==0)=0;
% 
% % Find the combined power loss associated with the main pump while also
% % accounting for the fact that since the main pump flow rate is different
% % than the HECM. This finds an "equivelant" power loss. If confused, think
% % about it from an energy perspective and relate the flows.
% 
% % The next lines assume...
% % When uncommented: That none of the flow back into a rail prevents losses.
% % When commented: That all of the flow back into a rail prevents losses.
% 
% % ************************************************************************
% % Q1_A_MPLoss(Q1_A_MPLoss<0) = 0;
% % Q1_B_MPLoss(Q1_B_MPLoss<0) = 0;
% % Q2_A_MPLoss(Q2_A_MPLoss<0) = 0;
% % Q2_B_MPLoss(Q2_B_MPLoss<0) = 0;
% % ************************************************************************
% 
% % % PLossMP1 = (abs(Q1_A.*PLossMP_A1)+abs(Q1_B.*PLossMP_B1))/(wMP*DMP);
% % % PLossMP2 = (abs(Q2_A.*PLossMP_A2)+abs(Q2_B.*PLossMP_B2))/(wMP*DMP);
% % % PLossMP3 = (abs(Q3_A.*PLossMP_A3)+abs(Q3_B.*PLossMP_B3))/(wMP*DMP);
% % % PLossMP4 = (abs(Q4_A.*PLossMP_A4)+abs(Q4_B.*PLossMP_B4))/(wMP*DMP);
% 
% PLossMP1 = ((Q1_A.*PLossMP_A1)+(Q1_B.*PLossMP_B1))/(wMP*DMP);
% PLossMP2 = ((Q2_A.*PLossMP_A2)+(Q2_B.*PLossMP_B2))/(wMP*DMP);
% PLossMP3 = ((Q3_A.*PLossMP_A3)+(Q3_B.*PLossMP_B3))/(wMP*DMP);
% PLossMP4 = ((Q4_A.*PLossMP_A4)+(Q4_B.*PLossMP_B4))/(wMP*DMP);
% 
% PLossMP = PLossMP1+PLossMP2;%+PLossMP3+PLossMP4;

MainPumpEfficiency = .885;

PLossMP = (1/MainPumpEfficiency-1)*(Q1_A.*CPR_A1+Q1_B.*CPR_B1+Q2_A.*CPR_A2+Q2_B.*CPR_B2+Q3_A.*CPR_A3+Q3_B.*CPR_B3+Q4_A.*CPR_A4+Q4_B.*CPR_B4);

%% Find "Constraint" Losses

% Find "Constraint Power Loss." Using the langrange multiplier method. The
% only constraint is that the net energy in the battery must be 0 over the
% drive cycle so that "free electrical energy" decisions are not made.

% Note that the energy generated and put into the battery will be less than
% the mechanical energy across the HECM due to losses.
% Find Electrical Power *GENERATED* when HECM is acting as a GENERATOR.
PLossConstraint_Generator1 = w1_Generator.*T1_Act*ElectricEff1;
PLossConstraint_Generator2 = w2_Generator.*T2_Act*ElectricEff2;
PLossConstraint_Generator3 = w3_Generator.*T3_Act*ElectricEff3;
PLossConstraint_Generator4 = w4_Generator.*T4_Act*ElectricEff4;

% Note that the energy consumed and taken from the battery will be more
% than the mechanical energy across the HECM due to losses.
% Find Electrical Power *CONSUMED* when HECM is acting as a MOTOR.
PLossConstraint_Motor1 = w1_Motor.*T1_Act/ElectricEff1;
PLossConstraint_Motor2 = w2_Motor.*T2_Act/ElectricEff2;
PLossConstraint_Motor3 = w3_Motor.*T3_Act/ElectricEff3;
PLossConstraint_Motor4 = w4_Motor.*T4_Act/ElectricEff4;

% PLossConstraints =  (PLossConstraint_Generator1+PLossConstraint_Generator2+PLossConstraint_Generator3+PLossConstraint_Generator4+PLossConstraint_Motor1+PLossConstraint_Motor2+PLossConstraint_Motor3+PLossConstraint_Motor4);
PLossConstraints =  PLossConstraint_Generator1+PLossConstraint_Motor1+PLossConstraint_Generator2+PLossConstraint_Motor2+PLossConstraint_Generator3+PLossConstraint_Motor3+PLossConstraint_Generator4+PLossConstraint_Motor4;

%% Sum All Losses/Costs ***************************************************

TotalLosses = PLossMP+PLoss_Hydraulic_HECM+PLoss_Electric_HECM;
Cost = PLossMP+PLoss_Hydraulic_HECM+PLoss_Electric_HECM+PLossConstraints*lambda1;

% Set losses for combinations that produce cavitation to inf so that they
% are not picked!
TotalLosses(P1_Cyl<0)=inf;
Cost(P1_Cyl<0)=inf;
TotalLosses(P2_Cyl<0)=inf;
Cost(P2_Cyl<0)=inf;
TotalLosses(P3_Cyl<0)=inf;
Cost(P3_Cyl<0)=inf;

%% Find Minimum

% Find the combination with the least power loss. "I" represents which
% combination produces the minimum power loss and M is the power loss for
% that combination.

[M, I] = min(Cost,[],2);

TotalCost = sum(M)*t_step*spacing;

toc

%% Check Constraint

% [DV_Ind_HECM1_Selected, PRA1_Selected, PRB1_Selected, DV_Ind_HECM2_Selected, PRA2_Selected, PRB2_Selected] = ind2sub(size(IndexMatrix),I);

% RealizedTotalLosses = zeros(1,length(t));
ConstraintCheck = zeros(1,length(t));

for i=1:length(t)

%     RealizedTotalLosses(i) = TotalLosses(i,I(i));
    ConstraintCheck(i) = PLossConstraints(i,I(i));

end

% Find the magnitude of the net energy of the battery divided by the sum 
% of the absolute value of the power acting on/from the system. (Then
% multiplied by 100 to get a percent.)

ConstraintCheckValue = 100*abs(sum(ConstraintCheck)/sum(abs(V1.*F1)))

BatteryCharge = -cumsum(ConstraintCheck*t_step);

figure
plot(t,BatteryCharge)
title('Battery Charge')
xlabel('Time [s])')
ylabel('Charge [J]')

BatterySize = max(BatteryCharge)-min(BatteryCharge)

%% %% %% %% %% %% %% Analysis %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%

% Section 1: Determine Selected Operation and Associated Metrics
% Section 2: CPR Flow Analysis
% Section 3: Energy Analysis and Efficiency
% Section 4: HECM Operation
% Section 5: 

%% Section 1: Determine Selected Operation and Associated Metrics

TotalLosses_Sel = zeros(1,length(t));

CPR_A1_Sel = zeros(1,length(t));
CPR_B1_Sel = zeros(1,length(t)); 
CPR_A2_Sel = zeros(1,length(t));
CPR_B2_Sel = zeros(1,length(t)); 
CPR_A3_Sel = zeros(1,length(t));
CPR_B3_Sel = zeros(1,length(t));
CPR_A4_Sel = zeros(1,length(t));
CPR_B4_Sel = zeros(1,length(t));

Q1_A_Sel = zeros(1,length(t));
Q1_B_Sel = zeros(1,length(t));
Q2_A_Sel = zeros(1,length(t));
Q2_B_Sel = zeros(1,length(t));
Q3_A_Sel = zeros(1,length(t));
Q3_B_Sel = zeros(1,length(t));
Q4_A_Sel = zeros(1,length(t));
Q4_B_Sel = zeros(1,length(t));

DeltaP_HECM1_Sel = zeros(1,length(t));
DeltaP_HECM2_Sel = zeros(1,length(t));
DeltaP_HECM3_Sel = zeros(1,length(t));
DeltaP_HECM4_Sel = zeros(1,length(t));

w1_Sel = zeros(1,length(t));
w2_Sel = zeros(1,length(t));
w3_Sel = zeros(1,length(t));
w4_Sel = W1;

CPR_F1_Sel = zeros(1,length(t));
CPR_F2_Sel = zeros(1,length(t));
CPR_F3_Sel = zeros(1,length(t));
CPR_T1_Sel = zeros(1,length(t));


E_HECM_Losses_Sel = zeros(1,length(t));
E_MP_Losses_Sel = zeros(1,length(t));

    T1_Act_Sel = zeros(length(t),1);
    T2_Act_Sel = zeros(length(t),1);
    T3_Act_Sel = zeros(length(t),1);
    T4_Act_Sel = zeros(length(t),1);
     
    Power1_Sel = zeros(length(t),1);
    Power2_Sel = zeros(length(t),1);
    Power3_Sel = zeros(length(t),1);
    Power4_Sel = zeros(length(t),1);

SelectedRailFlowMatrix = zeros(length(t),length(PR));

for i = 1:length(t)
    
    PLoss_Hydraulic_Sel = PLoss_Hydraulic_HECM(i,I(i));
    TotalLosses_Sel(i) = TotalLosses(i,I(i));
    E_HECM_Losses_Sel(i) = PLoss_Electric_HECM(i,I(i))+PLoss_Hydraulic_HECM(i,I(i));
    E_MP_Losses_Sel(i) = PLossMP(i,I(i));
    
    T1_Act_Sel(i) = T1_Act(i,I(i));
    T2_Act_Sel(i) = T2_Act(i,I(i));
    T3_Act_Sel(i) = T3_Act(i,I(i));
    T4_Act_Sel(i) = T4_Act(i,I(i));
    
    CPR_A1_Sel(i) = CPR_A1(i,I(i));
    CPR_B1_Sel(i) = CPR_B1(i,I(i));  
    CPR_A2_Sel(i) = CPR_A2(i,I(i));
    CPR_B2_Sel(i) = CPR_B2(i,I(i)); 
    CPR_A3_Sel(i) = CPR_A3(i,I(i));
    CPR_B3_Sel(i) = CPR_B3(i,I(i));
    CPR_A4_Sel(i) = CPR_A4(i,I(i));
    CPR_B4_Sel(i) = CPR_B4(i,I(i));
    
    Q1_A_Sel(i) = Q1_A(i,(I(i)));
    Q1_B_Sel(i) = Q1_B(i,(I(i)));
    Q2_A_Sel(i) = Q2_A(i,(I(i)));
    Q2_B_Sel(i) = Q2_B(i,(I(i)));
    Q3_A_Sel(i) = Q3_A(i,(I(i)));
    Q3_B_Sel(i) = Q3_B(i,(I(i)));
    Q4_A_Sel(i) = Q4_A(i,(I(i)));
    Q4_B_Sel(i) = Q4_B(i,(I(i)));
    
    DeltaP_HECM1_Sel(i) = DeltaP_HECM1(i,I(i));
    DeltaP_HECM2_Sel(i) = DeltaP_HECM2(i,I(i));
    DeltaP_HECM3_Sel(i) = DeltaP_HECM3(i,I(i));
    DeltaP_HECM4_Sel(i) = DeltaP_HECM4(i,I(i));
    
    w1_Sel(i) = w1(i,I(i));
    w2_Sel(i) = w2(i,I(i));
    w3_Sel(i) = w3(i,I(i));
    w4_Sel = W1;
    
    Power1_Sel(i) = w1_Sel(i)*T1_Act_Sel(i);
    Power2_Sel(i) = w2_Sel(i)*T2_Act_Sel(i);
    Power3_Sel(i) = w3_Sel(i)*T3_Act_Sel(i);
    Power4_Sel(i) = w4_Sel(i)*T4_Act_Sel(i);
    
    CPR_F1_Sel(i) = CPR_F1(I(i));
    CPR_F2_Sel(i) = CPR_F2(I(i));
    CPR_F3_Sel(i) = CPR_F3(I(i));
%     CPR_T1_Sel(i) = (CPR_A4_Sel(i)-CPR_B4_Sel(i))*D4*ScaleHECM4;
    CPR_T1_Sel(i) = CPR_T1(I(i));
    
    if CPR_A1_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q1_A_Sel(i);
    elseif CPR_A1_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q1_A_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q1_A_Sel(i);
    end
    
    if CPR_B1_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q1_B_Sel(i);
    elseif CPR_B1_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q1_B_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q1_B_Sel(i);
     end
    
    if CPR_A2_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q2_A_Sel(i);
    elseif CPR_A2_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q2_A_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q2_A_Sel(i);
    end
    
    if CPR_B2_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q2_B_Sel(i);
    elseif CPR_B2_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q2_B_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q2_B_Sel(i);
    end
    
    if CPR_A3_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q3_A_Sel(i);
    elseif CPR_A3_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q3_A_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q3_A_Sel(i);
    end
    
    if CPR_B3_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q3_B_Sel(i);
    elseif CPR_B3_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q3_B_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q3_B_Sel(i);
    end
    
    if CPR_A4_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q4_A_Sel(i);
    elseif CPR_A4_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q4_A_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q4_A_Sel(i);
    end
    
    if CPR_B4_Sel(i) == PR(1)
        SelectedRailFlowMatrix(i,1) = SelectedRailFlowMatrix(i,1)+Q4_B_Sel(i);
    elseif CPR_B4_Sel(i) == PR(2)
        SelectedRailFlowMatrix(i,2) = SelectedRailFlowMatrix(i,2)+Q4_B_Sel(i);
    else
        SelectedRailFlowMatrix(i,3) = SelectedRailFlowMatrix(i,3)+Q4_B_Sel(i);
    end
    
end

T1MAX = max(abs(T1_Act_Sel))
T2MAX = max(abs(T2_Act_Sel))
T3MAX = max(abs(T3_Act_Sel))
T4MAX = max(abs(T4_Act_Sel))

    Power1_MAX = max(abs(Power1_Sel))
    Power2_MAX = max(abs(Power2_Sel))
    [Power3_MAX PMax2_Ind] = max(abs(Power3_Sel))
    Power4_MAX = max(abs(Power4_Sel))
    
%% Section 2: Rail Flow Analysis

NetFlow_IntoTank = -cumsum(SelectedRailFlowMatrix(:,1)*t_step);
NetFlow_IntoR1 = -cumsum(SelectedRailFlowMatrix(:,2)*t_step);
% NetFlow_IntoR2 = -cumsum(SelectedRailFlowMatrix(:,3)*t_step);

figure
hold on
plot(t,NetFlow_IntoTank,t,NetFlow_IntoR1);%,t,NetFlow_IntoR2)
legend('Tank','CPR1')%,'CPR2')
title('Flow into CPRs')
ylabel('State of Accumulator [m^3]')
xlabel('Time [s]')

E_HECM_Losses_Total = sum(E_HECM_Losses_Sel)*t_step;
E_MP_Losses_Total = sum(E_MP_Losses_Sel)*t_step;

% if NetFlow_IntoR1(end)>0
%     E_MP_Losses_Total = E_MP_Losses_Total+NetFlow_IntoR1(end)*PR(2)*.11;
% end
% 
% if NetFlow_IntoR2(end)>0
%     E_MP_Losses_Total = E_MP_Losses_Total+2*NetFlow_IntoR2(end)*PR(3)*.12;
% end

%% Section 3: Energy Analysis and Efficiency

PositiveWorkPerformed1 = V1.*F1*t_step;
PositiveWorkPerformed1(PositiveWorkPerformed1>0)=0;
PositiveWorkPerformed2 = V2.*F2*t_step;
PositiveWorkPerformed2(PositiveWorkPerformed2>0)=0;
PositiveWorkPerformed3 = V3.*F3*t_step;
PositiveWorkPerformed3(PositiveWorkPerformed3>0)=0;
PositiveWorkPerformed4 = W1.*T1*t_step;
PositiveWorkPerformed4(PositiveWorkPerformed4>0)=0;
PositiveWorkPerformed = abs(sum(PositiveWorkPerformed1+PositiveWorkPerformed2+PositiveWorkPerformed3+PositiveWorkPerformed4));
E_Output_Positive = PositiveWorkPerformed;

RegenPotential_1 = V1.*F1*t_step;
RegenPotential_1(RegenPotential_1<0)=0;
RegenPotential_2 = V2.*F2*t_step;
RegenPotential_2(RegenPotential_2<0)=0;
RegenPotential_3 = V3.*F3*t_step;
RegenPotential_3(RegenPotential_3<0)=0;
RegenPotential_4 = T1.*W1*t_step;
RegenPotential_4(RegenPotential_4<0)=0;
RegenPotential = abs(sum(RegenPotential_1+RegenPotential_2+RegenPotential_3+RegenPotential_4));
E_Output_Negative = RegenPotential;

% E_Input_Net_time_BestCase = (CPR_A1_Sel.*Q1_A_Sel + CPR_B1_Sel.*Q1_B_Sel + CPR_A2_Sel.*Q2_A_Sel + CPR_B2_Sel.*Q2_B_Sel)*t_step + PLossMP_Sel*t_step;
% E_Input_Net_time_BestCase = (CPR_A1_Selected.*Q1_A_Selected + CPR_B1_Selected.*Q1_B_Selected + CPR_A2_Selected.*Q2_A_Selected + CPR_B2_Selected.*Q2_B_Selected)*t_step + TotalLosses_Selected*t_step;
% E_Input_Net_BestCase = sum(E_Input_Net_time_BestCase)
E_Input_Net_BestCase = sum(TotalLosses_Sel'-V1.*F1-V2.*F2-V3.*F3-T1.*W1)*t_step;%+NetFlow_IntoR1(end)*PR(2)*.11
E_Input_Net_BestCase2 = E_Output_Positive-E_Output_Negative+E_HECM_Losses_Total+E_MP_Losses_Total;

Efficiency_BestCase = 100*E_Output_Positive/E_Input_Net_BestCase
Efficiency_BestCase2 = 100*E_Output_Positive/E_Input_Net_BestCase2

return

%% Section 4: HECM Operation

PumpMotorEfficiencyMapping
hold on
plot(w1_Sel,DeltaP_HECM1_Sel,'bx','MarkerSize',1)
grid on
title('HECM 1 Operating Points 107cc Map')
xlabel('Angular Velocity [rad/s]')
ylabel('Pressure [Pa]')
hold off
%%
PumpMotorEfficiencyMapping
hold on
plot(w2_Sel,DeltaP_HECM2_Sel,'bx','MarkerSize',1)
grid on
title('HECM 2 Operating Points 107cc Map')
xlabel('Angular Velocity [rad/s]')
ylabel('Pressure [Pa]')
hold off
%%
PumpMotorEfficiencyMapping
hold on
plot(w3_Sel,DeltaP_HECM3_Sel,'bx','MarkerSize',1)
grid on
title('HECM 3 Operating Points 107cc Map')
xlabel('Angular Velocity [rad/s]')
ylabel('Pressure [Pa]')
hold off
%%
PumpMotorEfficiencyMapping
hold on
plot(w4_Sel,DeltaP_HECM4_Sel,'bx','MarkerSize',1)
grid on
title('HECM 4 Operating Points 107cc Map')
xlabel('Angular Velocity [rad/s]')
ylabel('Pressure [Pa]')
hold off

%% Section 5: HECM Selection Force Comparrisons (Linear Actuators)

fontsize = 10;

figure
plot(t,-F1)
hold on
plot(t,CPR_F1_Sel)
title('Actuator 1 Forces','FontSize',fontsize,'FontWeight','Normal')
% xlabel('Time [sec]')
% ylabel('Force [N]')
xlabel('Time [sec]','FontSize',fontsize,'FontWeight','Normal')
ylabel('Force [N]','FontSize',fontsize,'FontWeight','Normal')
% legend('Pressure Rail Force','Drive Cycle Force')
plot([zeros(1,length(CPR_F1)); t(end)*ones(1,length(CPR_F1))],[CPR_F1'; CPR_F1'],'k--')
set(gca,'FontSize',fontsize)
%%
figure
plot(t,-F2)
hold on
plot(t,CPR_F2_Sel)
title('Actuator 2 Forces','FontSize',fontsize,'FontWeight','Normal')
% xlabel('Time [sec]')
% ylabel('Force [N]')
xlabel('Time [sec]','FontSize',fontsize,'FontWeight','Normal')
ylabel('Force [N]','FontSize',fontsize,'FontWeight','Normal')
% legend('Pressure Rail Force','Drive Cycle Force')
plot([zeros(1,length(CPR_F2)); t(end)*ones(1,length(CPR_F2))],[CPR_F2'; CPR_F2'],'k')
set(gca,'FontSize',fontsize)
%%
figure
plot(t,-F3)
hold on
plot(t,CPR_F3_Sel)
title('Actuator 3 Forces','FontSize',fontsize,'FontWeight','Normal')
% xlabel('Time [sec]')
% ylabel('Force [N]')
xlabel('Time [sec]','FontSize',fontsize,'FontWeight','Normal')
ylabel('Force [N]','FontSize',fontsize,'FontWeight','Normal')
% legend('Pressure Rail Force','Drive Cycle Force')
plot([zeros(1,length(CPR_F3)); t(end)*ones(1,length(CPR_F3))],[CPR_F3'; CPR_F3'],'k')
set(gca,'FontSize',fontsize)
%%
figure
plot(t,-T1)
hold on
plot(t,CPR_T1_Sel)
title('Actuator 4 Torques','FontSize',fontsize,'FontWeight','Normal')
% xlabel('Time [sec]')
% ylabel('Force [N]')
xlabel('Time [sec]','FontSize',fontsize,'FontWeight','Normal')
ylabel('Torque [Nm]','FontSize',fontsize,'FontWeight','Normal')
% legend('Pressure Rail Torque','Drive Cycle Torque')
plot([zeros(1,length(CPR_T1)); t(end)*ones(1,length(CPR_T1))],[CPR_T1'; CPR_T1'],'k')
set(gca,'FontSize',fontsize)

%% Pie Chart

figure

subplot(1,2,1)
pie_Energy = [E_MP_Losses_Total E_HECM_Losses_Total E_Output_Positive];
pie(pie_Energy)
p = pie(pie_Energy);
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Main Pump Losses: ','HECM Losses: ','Output Work: '}; 
combinedtxt = strcat(txt,percentValues'); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
t = p(2);
t.FontSize = 14;
t = p(4);
t.FontSize = 14;
t = p(6);
t.FontSize = 14;
title('HHEA Losses','FontSize',16)

subplot(1,2,2)
pie_Energy = [RegenPotential E_Input_Net_BestCase2];
pie(pie_Energy)
q = pie(pie_Energy);
pText = findobj(q,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Regenerative Energy: ','"Input" Energy: '}; 
combinedtxt = strcat(txt,percentValues'); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
t = q(2);
t.FontSize = 14;
t = q(4);
t.FontSize = 14;
title('HHEA Energy Sources','FontSize',16)