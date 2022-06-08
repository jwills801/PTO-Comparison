% Run optimization with battery contraint
load JCB_LongerDC.mat
% masterScale = 1.25;
masterScale = 1;
% MaxT1_Act = inf;
% MaxT2_Act = inf;
% MaxT3_Act = inf;
% MaxT4_Act = inf;

maxRPM = 5000; %RPM %used to size the HECM

% (2 CPRS) - for NP
% MaxT1_Act = 24; 
% MaxT2_Act = 65; 
% MaxT3_Act = 36; 
% MaxT4_Act = 43; %torque limitation for each of the HECM modules

% (2 CPRS) - for C0
% MaxT1_Act = 24; 
% MaxT2_Act = 78; 
% MaxT3_Act = 39; 
% MaxT4_Act = 43; %torque limitation for each of the HECM modules

% LOWEST Possible 2 CPRs
% boom 24
% arm 65
% bkt 36
% sw 43

% (3 CPRs) - for NPP
% MaxT1_Act =13 %11 %12; 
% MaxT2_Act =22 %20 %23;
% MaxT3_Act =15 %15 %16;
% MaxT4_Act =22 %20 %22;


% (3 CPRs) - for C00
MaxT1_Act = 28; 
MaxT2_Act = 28;
MaxT3_Act = 28;
MaxT4_Act = 22;

% LOWEST Possible 3 CPRs
% boom 12
% arm 22
% bkt 16
% sw 22

%type of drive cycle operation:
% cycle = NinetyDeg;
cycle = Grading;
% cycle = Trenching;

load 28point1ccnolock_A6VM_4Q_20191210.mat

% Define an array of pressures
% PR = [0 1];
% PR_l = 0:0.01:0.1;
% PR_m = [0.51	0.505	0.5002	0.515	0.5296	0.525	0.5112	0.535	0.5216	0.5359	0.559];
% PR_h = ones(1,length(PR_l));
% for xx = 1:length(PR_l)
%     PR = [PR_l(xx) PR_m(xx) PR_h(xx)];
%     disp(PR)

% PR = [0.09 0.5359 1]; 
% PR = [0 0.5 1];
% PR = [0 1/3 2/3 1]; 
% PR = [0 1/4 1/2 3/4 1];

PR_A_decs = cell(4,2);
PR_B_decs = cell(4,2);
HECM_decs = cell(4,2);
Pcyl_set = cell(3,2);
% elim = cell(3,2);
HECM_torque = cell(4,2);
HECM_angvel = cell(4,2);

for p_rail = 1
    if p_rail == 1
        PR = [0 0.5 1]*29e6;
        all_options = {'C00'};
        
        MaxT1_Act =28 %11 %12; 
        MaxT2_Act =28 %20 %23;
        MaxT3_Act =28 %15 %16;
        MaxT4_Act =22 %20 %22;
    else
        PR = [0.09 0.5359 1]*32e6; 
        all_options = {'C0P'};
        
        MaxT1_Act =13 %11 %12; 
        MaxT2_Act =22 %20 %23;
        MaxT3_Act =15 %15 %16;
        MaxT4_Act =22 %20 %22;
    end
 

    Make_losses_general;
    %return
    %%
    % 2 rail option
    %     all_options = {'NP'}
    %all_options = {'NP','NM','C0'}


    % % 3 rail options
%     all_options = {'C0P'};
    % all_options = { 'NPP';
    %                 'C0P';
    %                 'NMP';
    %                 'CP0';
    %                 'C00';
    %                 'CM0';
    %                 'NPM';
    %                 'C0M';
    %                 'NMM'}

    % 4 rail options
    % all_options = {'NPPP','NPPM','NPMP', 'NMPP','NPMM','NMPM','NMMP',...
    %     'C0PP','CP0P','CPP0',...
    %     'C0PM','C0MP','CP0M','CM0P','CMP0','CPM0',...
    %     'C00P','C0P0','CP00'}
    % Each region
    %all_options = {'NPPP','C0PP','CP0P','CPP0','C00P*','C0P0','CP00'}
    %all_options = {'NMPP','C0PP','CM0P','CMP0','C00P*','C0P0'}
    %all_options = {'NMMP','C0MP','CM0P','C00P*'}
    %all_options = {'NPPM','C0PM','CP0M','CPP0','C0P0','CP00'} %=> 'no solution'
    %all_options = {'NPMP','C0MP','CP0P','CPM0','C00P*','CP00'} 
    %all_options = {'NPMM','CP0M','CPM0','CP00'} %=> 'no solution'
    %all_options = {'NMPM','C0PM','CMP0','C0P0'} %=> 'no solution'

    % 5 rails
    % all_options = {'NPPPP','C0PPP','CP0PP','CPP0P','CPPP0',...
    %    'C00PP','C0P0P','C0PP0','CP00P','CP0P0','CPP00',...
    %    'C000P','C00P0','C0P00','CP000'}

    %all_options = {'NPPPP'... 1
    %             'C0PPP'... 2
    %             'NMPPP'... 3
    %             'CP0PP'... 4
    %             'C00PP'... 5
    %             'CM0PP'... 6
    %             'NPMPP'... 7

    %             'C0MPP'... 8
    %             'NMMPP'... 9
    %             'CPP0P'... 10
    %             'C0P0P'... 11
    %             'CMP0P'... 12
    %             'CP00P'... 13
    %             'C000P'... 14
    %             'CM00P'... 15 *
    %             'CPM0P'... 16
    %             'C0M0P'... 17
    %             'CMM0P'... 18
    %             'NPPMP'... 19
    %             'C0PMP'... 20 
    %             'NMPMP'... 21
    %             'CP0MP'... 22
    %             'C00MP'... 23
    %             'CM0MP'... 24
    %             'NPMMP'... 25
    %             'C0MMP'... 26
    %             'NMMMP'... 27
    %             'CPPP0'... 28
    %             'C0PP0'... 29
    %             'CMPP0'... 30
    %             'CP0P0'... 31
    %             'C00P0'... 32
    %             'CM0P0'... 33
    %             'CPMP0'... 34
    %             'C0MP0'... 35
    %             'CMMP0'... 36
    %             'CPP00'... 37
    %             'C0P00'... 38
    %             'CMP00'... 39
    %             'CP000'... 40 *
    %             'CM000'... 41
    %             'CPM00'... 42
    %             'C0M00'... 43
    %             'CMM00'... 44
    %             'CPPM0'... 45
    %             'C0PM0'... 46
    %             'CMPM0'... 47
    %             'CP0M0'... 48
    %             'C00M0'... 49 *
    %             'CM0M0'... 50
    %             'CPMM0'... 51
    %             'C0MM0'... 52
    %             'CMMM0'... 53
    %             'NPPPM'... 54
    %             'C0PPM'... 55
    %             'NMPPM'... 56
    %             'CP0PM'... 57
    %             'C00PM'... 58
    %             'CM0PM'... 59
    %             'NPMPM'... 60
    %             'C0MPM'... 61
    %             'NMMPM'... 62
    %             'CPP0M'... 63
    %             'C0P0M'... 64
    %             'CMP0M'... 65
    %             'CP00M'... 66
    %             'C000M'... 67
    %             'CM00M'... 68 *
    %             'CPM0M'... 69
    %             'C0M0M'... 70
    %             'CMM0M'... 71
    %             'NPPMM'... 72
    %             'C0PMM'... 73
    %             'NMPMM'... 74
    %             'CP0MM'... 75
    %             'C00MM'... 76
    %             'CM0MM'... 77
    %             'NPMMM'... 78
    %             'C0MMM'... 79
    %             'NMMMM'}; %80


    % all_options = {'NPPPP'};
    % all_options = {'CM00P'};        

    tic
    Eff_all=[];
    NetRail_all=[];
    for case_no = 1:length(all_options)
    options = all_options{case_no}
    [MPstr,Jstr,BATstr,lam]=expand_options(options,PR)
    eval(MPstr);
    eval(BATstr);
    eval(['f=@(lam) -sum(',Jstr,');']);
    options = optimset('PlotFcns',@optimplotfval);
    tic
    lam = fminsearch(f,lam,options)
    toc
    %return 
    %% evaluate the optimal solution
    eval(['[cost,d]=',Jstr,';'])
    disp('Cost is: ')
    cost;
    sz = size(HECMLosses); 
    I_t = [1:sz(1)]';
    dt = mean(diff(t));
    d_ind = sub2ind(sz,I_t,d);

    for kk=1:length(PR)
    NetRailEnergy(kk)=sum(QR{kk}(d_ind))*PR(kk)*dt;
    end

    bat_act = battery_power_HECM;

    for kk=2:length(PR)
    if NetRailEnergy(kk)>0
        bat_act = bat_act + battery_power_rail_p{kk-1};
    elseif NetRailEnergy(kk)<0
        bat_act = bat_act + battery_power_rail_m{kk-1};
    end
    end

    bat_act = bat_act(d_ind);

    NetRailEnergy;
    ActualMPLoss = NetRailEnergy.*((NetRailEnergy>0).*[0,R_P_loss'] ...
                             - (NetRailEnergy<0).*[0,R_M_loss']);
    TotalMPWork = sum(NetRailEnergy+ActualMPLoss);
    Netbattery = sum(bat_act)*dt
    NetWorkout = -sum(F1.*V1+F2.*V2+F3.*V3+T4.*W4)*dt;
    Totallosses = sum(HECMLosses(d_ind))*dt+sum(ActualMPLoss);
    tmp=F1.*V1; PosWork=-sum(tmp(tmp<0))*dt;
    tmp=F2.*V2; PosWork=PosWork-sum(tmp(tmp<0))*dt;
    tmp=F3.*V3; PosWork=PosWork-sum(tmp(tmp<0))*dt;
    tmp=T4.*W4; PosWork=PosWork-sum(tmp(tmp<0))*dt;
    tmp=F1.*V1; NegWork=sum(tmp(tmp>0))*dt;
    tmp=F2.*V2; NegWork=NegWork+sum(tmp(tmp>0))*dt;
    tmp=F3.*V3; NegWork=NegWork+sum(tmp(tmp>0))*dt;
    tmp=T4.*W4; NegWork=NegWork+sum(tmp(tmp>0))*dt;
    Efficiency = PosWork/(NetWorkout+Totallosses)
    [TotalMPWork,NetWorkout+Totallosses-Netbattery]
    [sum(cost)*dt,Totallosses]

    NetRail_all(case_no,:) = NetRailEnergy;
    Eff_all(case_no)=Efficiency
    %     case_eff(xx) = Efficiency;

    figure(2);
    subplot(211); plot(t, cumsum(bat_act)*dt); c_battery = cumsum(bat_act)*dt;
    title(['Battery end state = ',num2str(c_battery(end)),'J'])
    if  length(PR)==5
    subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind),QR{5}(d_ind)])*dt); 
    legend('QR{1}','QR{2}','QR{3}','QR{4}','QR{5}');
    elseif length(PR)==4
    subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind)])*dt); 
    legend('QR{1}','QR{2}','QR{3}','QR{4}');
    elseif length(PR)==3
    subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind)])*dt); 
    legend('QR{1}','QR{2}','QR{3}');
    elseif length(PR)==2
    subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind)])*dt); 
    legend('QR{1}','QR{2}');
    end
    ylabel('Cumulative Rail Flow');
    xlabel('Time - s');
    %end
    %save SC_5rails_results.mat NetRail_all Eff_all PR ScaleMaxT optimal
    %toc
    %% Check force/torque

    %     figure(3);
    %     subplot(221); plot(t, -F1, '-', t([1,end]),ones(2,1)*Frange1', '-', t, CPR_F1(d),'-.','linewidth',2); ylabel('Force1 N'); xlabel('Time - s');grid
    %     subplot(222); plot(t, -F2, '-', t([1,end]),ones(2,1)*Frange2', '-', t, CPR_F2(d),'-.','linewidth',2); ylabel('Force2 N'); xlabel('Time - s');grid
    %     subplot(223); plot(t, -F3, '-', t([1,end]),ones(2,1)*Frange3', '-', t, CPR_F3(d),'-.','linewidth',2); ylabel('Force3 N'); xlabel('Time - s');grid
    %     subplot(224); plot(t, -T4, '-', t([1,end]),ones(2,1)*Trange4', '-', t, CPR_T1(d),'-.','linewidth',2); ylabel('Torque4 Nm'); xlabel('Time - s');grid

    %% Check flow/speed
    %     figure(4);
    %     subplot(221); plot(t, V1, '-', t, -Q_HECM_1(d_ind)/ARod1,'-.'); ylabel('Speed1 m/s'); xlabel('Time - s');grid
    %     subplot(222); plot(t, V2, '-', t, -Q_HECM_2(d_ind)/ARod2,'-.'); ylabel('Speed2 m/s'); xlabel('Time - s');grid
    %     subplot(223); plot(t, V3, '-', t, -Q_HECM_3(d_ind)/ARod3,'-.'); ylabel('Speed3 m/s'); xlabel('Time - s');grid
    %     subplot(224); plot(t, W4, '-', t, w4(d_ind),'-.'); ylabel('Speed4 rad/s'); xlabel('Time - s');grid

    %% Check HECM E-motor operating points - (w,T)
    %     figure(5)
    %     subplot(211); plot(w1(d_ind), T1_Act(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]')
    %     subplot(212); plot(w2(d_ind), T2_Act(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]')
    max_power_peak_cornor = ...
    [max(abs(w1(d_ind).*T1_Act(d_ind))),max(abs(w1(d_ind))*max(abs(T1_Act(d_ind))));
     max(abs(w2(d_ind).*T2_Act(d_ind))),max(abs(w2(d_ind))*max(abs(T2_Act(d_ind))));
     max(abs(w3(d_ind).*T3_Act(d_ind))),max(abs(w3(d_ind))*max(abs(T3_Act(d_ind))));
     max(abs(w4(d_ind).*T4_Act(d_ind))),max(abs(w4(d_ind))*max(abs(T4_Act(d_ind))));];

    %% Check HECM P/M operating points - (w,P) vs HECM PM efficiency
    %Mech=T1_Act_Mapping.*w1rad_Mapping;
    %Hyd=Q1_Act_Mapping.*P1_Mapping;
    %Eff=Hyd./(Mech+eps).*(Hyd>0)+Mech./(Hyd+eps).*(Hyd<0);
    disp('Calculating HECM hydraulic efficiencies')
    HydHECM=DeltaP_HECM1(d_ind).*Q_HECM_1(d_ind);
    MechHECM=w1(d_ind).*T1_Act(d_ind);    
    indP=HydHECM>0; eff(1,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
%     indP = MechHECM>0; eff(1,1) = sum(HydHECM(indP))/sum(MechHECM(indP));
    xx1_P = HydHECM(indP);
    xx2_P = MechHECM(indP);
    indM=HydHECM<0; eff(1,2)=sum(MechHECM(indM))/sum(HydHECM(indM));
%     indM=MechHECM<0; eff(1,2)=sum(MechHECM(indM))/sum(HydHECM(indM));
    xx1_M = HydHECM(indM);
    xx2_M = MechHECM(indM);
    HydHECM=DeltaP_HECM2(d_ind).*Q_HECM_2(d_ind);
    MechHECM=w2(d_ind).*T2_Act(d_ind);
    %indP=HydHECM>0; eff(2,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
    indP=MechHECM>0; eff(2,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
    %indM=HydHECM<0; eff(2,2)=sum(MechHECM(indM))/sum(HydHECM(indM));
    indM=MechHECM<0; eff(2,2)=sum(MechHECM(indM))/sum(HydHECM(indM));
    
    a1 = zeros(size(T1_Act));
    a1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0) = DeltaP_HECM1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0).*Q_HECM_1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0);
    a1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0) = -T1_Act(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0).*w1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0);

    b1 = T1_Act.*w1-DeltaP_HECM1.*Q_HECM_1;
    b1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0) = T1_Act(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0).*w1(T1_Act.*w1>0 & DeltaP_HECM1.*Q_HECM_1>0);
    b1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0) = -DeltaP_HECM1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0).*Q_HECM_1(T1_Act.*w1<0 & DeltaP_HECM1.*Q_HECM_1<0);

    a2 = zeros(size(T2_Act));
    a2(T2_Act.*w2>0 & DeltaP_HECM2.*Q_HECM_2>0) = DeltaP_HECM2(T2_Act.*w2>0 & DeltaP_HECM2.*Q_HECM_2>0).*Q_HECM_2(T2_Act.*w2>0 & DeltaP_HECM2.*Q_HECM_2>0);
    a2(T2_Act.*w2<0 & DeltaP_HECM2.*Q_HECM_2<0) = -T2_Act(T2_Act.*w2<0 & DeltaP_HECM2.*Q_HECM_2<0).*w2(T2_Act.*w2<0 & DeltaP_HECM2.*Q_HECM_2<0);

    b2 = T2_Act.*w2-DeltaP_HECM2.*Q_HECM_2;
    b2(T2_Act.*w2>0 & DeltaP_HECM2.*Q_HECM_2>0) = T2_Act(T2_Act.*w2>0 & DeltaP_HECM2.*Q_HECM_2>0).*w2(T2_Act.*w2>0 & DeltaP_HECM2.*Q_HECM_2>0);
    b2(T2_Act.*w2<0 & DeltaP_HECM2.*Q_HECM_2<0) = -DeltaP_HECM2(T2_Act.*w2<0 & DeltaP_HECM2.*Q_HECM_2<0).*Q_HECM_2(T2_Act.*w2<0 & DeltaP_HECM2.*Q_HECM_2<0);


    EnergyConversionHECM_Percent = [sum(a1(d_ind))/sum(b1(d_ind)), sum(a2(d_ind))/sum(b2(d_ind))]*100;


    %     if max(eff(:)>1)
    %         disp('HECM Hydraulic Efficiency > 100%')
    %         return
    %     end

        figure(p_rail + 6); clf;
        subplot(211); plot(w1(d_ind), DeltaP_HECM1(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]');
        hold on; [c1,h1]=contour(w1rad_Mapping,P1_Mapping,Eff,[-0.1:0.1:1]); clabel(c1,h1)
        axis(1.05*[min(w1(d_ind)),max(w1(d_ind)),min(DeltaP_HECM1(d_ind)),max(DeltaP_HECM1(d_ind))])
        title(['Pumping/Motoring efficiency: ',num2str(eff(1,:)*100),'%'])
    
        subplot(212); plot(w2(d_ind), DeltaP_HECM2(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]')
        hold on; [c2,h2]=contour(w1rad_Mapping,P1_Mapping,Eff,[-0.1:0.1:1]); clabel(c2,h2)
        axis(1.05*[min(w2(d_ind)),max(w2(d_ind)),min(DeltaP_HECM2(d_ind)),max(DeltaP_HECM2(d_ind))])
        title(['Pumping/Motoring efficiency: ',num2str(eff(2,:)*100),'%'])
    %     [TotalMPWork,NetWorkout+Totallosses]
    EnergyConservation = [Netbattery+NetRail_all(1)-(Totallosses+NetWorkout)]
    end
    NetRail_all
    Eff_list = Eff_all'

    [PosWork/1000 sum(ActualMPLoss)/1000 sum(HECMLosses(d_ind))*dt/1000]
    [NegWork/1000 (NetWorkout+Totallosses)/1000]

    %% Analysing the decisions made:
    PRA1_mod = repmat(PRA1',size(ones_t));
    PRB1_mod = repmat(PRB1',size(ones_t));

    PRA2_mod = repmat(PRA2',size(ones_t));
    PRB2_mod = repmat(PRB2',size(ones_t));

    PRA3_mod = repmat(PRA3',size(ones_t));
    PRB3_mod = repmat(PRB3',size(ones_t));

    PRA4_mod = repmat(PRA4',size(ones_t));
    PRB4_mod = repmat(PRB4',size(ones_t));

    % finding the pressure rail selections and adding them to a cell array:
    
    PRA1_mod_dec = PRA1_mod(d_ind);
    PR_A_decs{1,p_rail} = PRA1_mod_dec;
    
    PRB1_mod_dec = PRB1_mod(d_ind);
    PR_B_decs{1,p_rail} = PRB1_mod_dec;
    
    PRA2_mod_dec = PRA2_mod(d_ind);
    PR_A_decs{2,p_rail} = PRA2_mod_dec;
    
    PRB2_mod_dec = PRB2_mod(d_ind);
    PR_B_decs{2,p_rail} = PRB2_mod_dec;
    
    PRA3_mod_dec = PRA3_mod(d_ind);
    PR_A_decs{3,p_rail} = PRA3_mod_dec;
    
    PRB3_mod_dec = PRB3_mod(d_ind);
    PR_B_decs{3,p_rail} = PRB3_mod_dec;
    
    PRA4_mod_dec = PRA4_mod(d_ind);
    PR_A_decs{4,p_rail} = PRA4_mod_dec;
    
    PRB4_mod_dec = PRB4_mod(d_ind);
    PR_B_decs{4,p_rail} = PRB4_mod_dec;
    
    % HECM Decisions:
    %delta pressures - HECM
    HECM_decs{1,p_rail} = DeltaP_HECM1(d_ind);
    HECM_decs{2,p_rail} = DeltaP_HECM2(d_ind);
    HECM_decs{3,p_rail} = DeltaP_HECM3(d_ind);
    HECM_decs{4,p_rail} = DeltaP_HECM4(d_ind);
    
    %torques - HECM
    HECM_torque{1,p_rail} = T1_Act(d_ind);
    HECM_torque{2,p_rail} = T2_Act(d_ind);
    HECM_torque{3,p_rail} = T3_Act(d_ind);
    HECM_torque{4,p_rail} = T4_Act(d_ind);
    
    %angular velocities - HECM
    HECM_angvel{1,p_rail} = w1(d_ind);
    HECM_angvel{2,p_rail} = w2(d_ind);
    HECM_angvel{3,p_rail} = w3(d_ind);
    HECM_angvel{4,p_rail} = w4(d_ind);
    
    %storing the cylinder pressures:
    Pcyl_set{1,p_rail} = P1_Cyl;
    Pcyl_set{2,p_rail} = P2_Cyl;
    Pcyl_set{3,p_rail} = P3_Cyl;
    %% Plotting the bar graphs:
    X = categorical({'PR(1)','PR(2)','PR(3)'});
    X = reordercats(X,{'PR(1)','PR(2)','PR(3)'});
%     figure(8)
% 
%     subplot(2,2,1);
%     bar(X,[length(find(PR(1) == PRB1_mod_dec)),length(find(PR(2) == PRB1_mod_dec)),length(find(PR(3) == PRB1_mod_dec))],0.3)
%     subplot(2,2,2);
%     bar(X,[length(find(PR(1) == PRB2_mod_dec)),length(find(PR(2) == PRB2_mod_dec)),length(find(PR(3) == PRB2_mod_dec))],0.3)
%     subplot(2,2,3);
%     bar(X,[length(find(PR(1) == PRB3_mod_dec)),length(find(PR(2) == PRB3_mod_dec)),length(find(PR(3) == PRB3_mod_dec))],0.3)
%     subplot(2,2,4);
%     bar(X,[length(find(PR(1) == PRB4_mod_dec)),length(find(PR(2) == PRB4_mod_dec)),length(find(PR(3) == PRB4_mod_dec))],0.3)

    %EnergyConservation = [Netbattery-(Totallosses+NetWorkout)]

    % sum([PosWork/1000 sum(ActualMPLoss)/1000 sum(HECMLosses(d_ind))*dt/1000]) / sum([NegWork/1000 NetWorkout+Totallosses/1000])
    % end
    decisions(:,p_rail) = d_ind;
%  return
end
return
% figure(9)
% plot(t_all,PR_A_decs{1}')
% hold on
% plot(t_all,PR_A_decs{5}')
% title('Comparison of pressure rail selections for Actuator one - Cap side')
% legend('Original PR','Modified PR')
% 
% figure(10)
% plot(t_all,PR_B_decs{1}')
% hold on
% plot(t_all,PR_B_decs{5}')
% title('Comparison of pressure rail selections for Actuator one - Rod side')
% legend('Original PR','Modified PR')
% 
% figure(11)
% plot(t_all,PR_A_decs{2}')
% hold on
% plot(t_all,PR_A_decs{6}')
% title('Comparison of pressure rail selections for Actuator two - Cap side')
% legend('Original PR','Modified PR')
% 
% figure(12)
% plot(t_all,PR_B_decs{2}')
% hold on
% plot(t_all,PR_B_decs{6}')
% title('Comparison of pressure rail selections for Actuator two - Rod side')
% legend('Original PR','Modified PR')

%% Rail Forces Comparison:

RF_org_1 = PR_A_decs{1}*ACap1 - PR_B_decs{1}*ARod1;
RF_mod_1 = PR_A_decs{5}*ACap1 - PR_B_decs{5}*ARod1;
% figure(13)
% plot([1:2000],RF_org_1([1:2000]))
% hold on
% plot([1:2000],RF_mod_1([1:2000]))
% hold on
% plot([1:2000], -F1([1:2000]))
% legend('Original PR','Modified PR','Actual Force')
% xlabel('Time')
% ylabel('Force (N)')
% %% HECM Forces Comparison:
% 
% %defining the net force for the two different cases:
% Net_F_org_1 = PR_A_decs{1}*ACap1 - (PR_B_decs{1}+HECM_decs{1})*ARod1;
% Net_F_mod_1 = PR_A_decs{5}*ACap1 - (PR_B_decs{5}+HECM_decs{5})*ARod1;
% 
% %finding out percentage contribution of HECM to the net force:
% HECM_org_pct_1 = (HECM_decs{1}*ARod1./Net_F_org_1)*100;
% HECM_mod_pct_1 = (HECM_decs{5}*ARod1./Net_F_mod_1)*100;
%     
% HECM_pct_diff_1 = (abs(HECM_org_pct_1) - abs(HECM_mod_pct_1));
% 
% %% Using cylinder pressure for analysis:
for mm = 1:3
    Pcyl_selected_1 = Pcyl_set{mm}; %selecting the set corresponding to the original PR
    Pcyl_selected_2 = Pcyl_set{mm + 3}; %selecting the set corresponding to the new PR
    for nn = 1:length(t_all)
        elim_1(nn) = length(find(Pcyl_selected_1(nn,:) < -1e5)); %finding the options eliminated bcoz of cavitation in og PR
        elim_2(nn) = length(find(Pcyl_selected_2(nn,:) < -1e5)); %finding the options eliminated bcoz of cavitation in new PR
        if mm == 1
            to_be_checked_1 = find(elim_1 ~= elim_2);
        elseif mm == 2
            to_be_checked_2 = find(elim_1 ~= elim_2);
        else
            to_be_checked_3 = find(elim_1 ~= elim_2);
        end
    end
end
% 
% %% Using the time instants where disproportionate options were eliminated:
% 
%for actuator 1:
key_dec_org_cap_1 = PR_A_decs{1}(to_be_checked_1');
key_dec_org_rod_1 = PR_B_decs{1}(to_be_checked_1');

key_dec_mod_cap_1 = PR_A_decs{5}(to_be_checked_1');
key_dec_mod_rod_1 = PR_B_decs{5}(to_be_checked_1');

key_dec_org_1 = [to_be_checked_1' key_dec_org_cap_1 key_dec_org_rod_1];
key_dec_mod_1 = [to_be_checked_1' key_dec_mod_cap_1 key_dec_mod_rod_1];

%for actuator 2:
key_dec_org_cap_2 = PR_A_decs{2}(to_be_checked_2');
key_dec_org_rod_2 = PR_B_decs{2}(to_be_checked_2');

key_dec_mod_cap_2 = PR_A_decs{6}(to_be_checked_2');
key_dec_mod_rod_2 = PR_B_decs{6}(to_be_checked_2');

key_dec_org_2 = [to_be_checked_2' key_dec_org_cap_2 key_dec_org_rod_2];
key_dec_mod_2 = [to_be_checked_2' key_dec_mod_cap_2 key_dec_mod_rod_2];

%for actuator 3:
key_dec_org_cap_3 = PR_A_decs{3}(to_be_checked_3');
key_dec_org_rod_3 = PR_B_decs{3}(to_be_checked_3');

key_dec_mod_cap_3 = PR_A_decs{7}(to_be_checked_3');
key_dec_mod_rod_3 = PR_B_decs{7}(to_be_checked_3');

key_dec_org_3 = [to_be_checked_3' key_dec_org_cap_3 key_dec_org_rod_3];
key_dec_mod_3 = [to_be_checked_3' key_dec_mod_cap_3 key_dec_mod_rod_3];
% 
% %% showing how the decisions change:
% 
for i = 1:length(t_all)
    %for actuator 1:
    
    %original pressure rail set cap side:
    if PR_A_decs{1}(i) == 0*Pmax
        PR_A_decs_simple{1}(i) = 1;
    elseif PR_A_decs{1}(i) == 0.5*Pmax
        PR_A_decs_simple{1}(i) = 2;
    elseif PR_A_decs{1}(i) == 1*Pmax
        PR_A_decs_simple{1}(i) = 3;
    end
    
    %original pressure rail set rod side:
    if PR_B_decs{1}(i) == 0*Pmax
        PR_B_decs_simple{1}(i) = 1;
    elseif PR_B_decs{1}(i) == 0.5*Pmax
        PR_B_decs_simple{1}(i) = 2;
    elseif PR_B_decs{1}(i) == 1*Pmax
        PR_B_decs_simple{1}(i) = 3;
    end
    
    %modified pressure rail set cap side:
    if PR_A_decs{5}(i) == 0.09*Pmax
        PR_A_decs_simple{5}(i) = 1;
    elseif PR_A_decs{5}(i) == 0.5359*Pmax
        PR_A_decs_simple{5}(i) = 2;
    elseif PR_A_decs{5}(i) == 1*Pmax
        PR_A_decs_simple{5}(i) = 3;
    end
    
    %modified pressure rail set rod side:
    if PR_B_decs{5}(i) == 0.09*Pmax
        PR_B_decs_simple{5}(i) = 1;
    elseif PR_B_decs{5}(i) == 0.5359*Pmax
        PR_B_decs_simple{5}(i) = 2;
    elseif PR_B_decs{5}(i) == 1*Pmax
        PR_B_decs_simple{5}(i) = 3;
    end
    
    %for actuator 2:
    %original pressure rail set cap side:
    if PR_A_decs{2}(i) == 0*Pmax
        PR_A_decs_simple{2}(i) = 1;
    elseif PR_A_decs{2}(i) == 0.5*Pmax
        PR_A_decs_simple{2}(i) = 2;
    elseif PR_A_decs{2}(i) == 1*Pmax
        PR_A_decs_simple{2}(i) = 3;
    end
    
    %original pressure rail set rod side:
    if PR_B_decs{2}(i) == 0*Pmax
        PR_B_decs_simple{2}(i) = 1;
    elseif PR_B_decs{2}(i) == 0.5*Pmax
        PR_B_decs_simple{2}(i) = 2;
    elseif PR_B_decs{2}(i) == 1*Pmax
        PR_B_decs_simple{2}(i) = 3;
    end
    
    %modified pressure rail set cap side:
    if PR_A_decs{6}(i) == 0.09*Pmax
        PR_A_decs_simple{6}(i) = 1;
    elseif PR_A_decs{6}(i) == 0.5359*Pmax
        PR_A_decs_simple{6}(i) = 2;
    elseif PR_A_decs{6}(i) == 1*Pmax
        PR_A_decs_simple{6}(i) = 3;
    end
    
    %modified pressure rail set rod side:
    if PR_B_decs{6}(i) == 0.09*Pmax
        PR_B_decs_simple{6}(i) = 1;
    elseif PR_B_decs{6}(i) == 0.5359*Pmax
        PR_B_decs_simple{6}(i) = 2;
    elseif PR_B_decs{6}(i) == 1*Pmax
        PR_B_decs_simple{6}(i) = 3;
    end
    
    %for actuator 3:
    %original pressure rail set cap side:
    if PR_A_decs{3}(i) == 0*Pmax
        PR_A_decs_simple{3}(i) = 1;
    elseif PR_A_decs{3}(i) == 0.5*Pmax
        PR_A_decs_simple{3}(i) = 2;
    elseif PR_A_decs{3}(i) == 1*Pmax
        PR_A_decs_simple{3}(i) = 3;
    end
    
    %original pressure rail set rod side:
    if PR_B_decs{3}(i) == 0*Pmax
        PR_B_decs_simple{3}(i) = 1;
    elseif PR_B_decs{3}(i) == 0.5*Pmax
        PR_B_decs_simple{3}(i) = 2;
    elseif PR_B_decs{3}(i) == 1*Pmax
        PR_B_decs_simple{3}(i) = 3;
    end
    
    %modified pressure rail set cap side:
    if PR_A_decs{7}(i) == 0.09*Pmax
        PR_A_decs_simple{7}(i) = 1;
    elseif PR_A_decs{7}(i) == 0.5359*Pmax
        PR_A_decs_simple{7}(i) = 2;
    elseif PR_A_decs{7}(i) == 1*Pmax
        PR_A_decs_simple{7}(i) = 3;
    end
    
    %modified pressure rail set rod side:
    if PR_B_decs{7}(i) == 0.09*Pmax
        PR_B_decs_simple{7}(i) = 1;
    elseif PR_B_decs{7}(i) == 0.5359*Pmax
        PR_B_decs_simple{7}(i) = 2;
    elseif PR_B_decs{7}(i) == 1*Pmax
        PR_B_decs_simple{7}(i) = 3;
    end
    
    %for actuator 4:
    %original pressure rail set cap side:
    if PR_A_decs{4}(i) == 0*Pmax
        PR_A_decs_simple{4}(i) = 1;
    elseif PR_A_decs{4}(i) == 0.5*Pmax
        PR_A_decs_simple{4}(i) = 2;
    elseif PR_A_decs{4}(i) == 1*Pmax
        PR_A_decs_simple{4}(i) = 3;
    end
    
    %original pressure rail set rod side:
    if PR_B_decs{4}(i) == 0*Pmax
        PR_B_decs_simple{4}(i) = 1;
    elseif PR_B_decs{4}(i) == 0.5*Pmax
        PR_B_decs_simple{4}(i) = 2;
    elseif PR_B_decs{4}(i) == 1*Pmax
        PR_B_decs_simple{4}(i) = 3;
    end
    
    %modified pressure rail set cap side:
    if PR_A_decs{8}(i) == 0.09*Pmax
        PR_A_decs_simple{8}(i) = 1;
    elseif PR_A_decs{8}(i) == 0.5359*Pmax
        PR_A_decs_simple{8}(i) = 2;
    elseif PR_A_decs{8}(i) == 1*Pmax
        PR_A_decs_simple{8}(i) = 3;
    end
    
    %modified pressure rail set rod side:
    if PR_B_decs{8}(i) == 0.09*Pmax
        PR_B_decs_simple{8}(i) = 1;
    elseif PR_B_decs{8}(i) == 0.5359*Pmax
        PR_B_decs_simple{8}(i) = 2;
    elseif PR_B_decs{8}(i) == 1*Pmax
        PR_B_decs_simple{8}(i) = 3;
    end
end
% figure(14)
% plot(t_all,PR_B_decs_simple{1},t_all,PR_B_decs_simple{5})
% ylim([0 4])
% set(gca,'ytick',0:1:4)
% xlabel('Time')
% ylabel('Rail Selected')
% legend('Old PRs','New PRs')
% title('Pressure Rail Selections for Actuator 1 - Rod Side')
% % 
% figure(15)
% subplot(2,1,1);
% plot([322:409],PR_B_decs_simple{1}([322:409]),[322:409],PR_B_decs_simple{5}([322:409]))
% ylim([0 4])
% set(gca,'ytick',0:1:4)
% xlabel('Time')
% ylabel('Rail Selected')
% legend('Old PRs','New PRs')
% title('Pressure Rail Selections for Actuator 1 - Rod Side')
% subplot(2,1,2);
% plot([322:409],PR_A_decs_simple{1}([322:409]),[322:409],PR_A_decs_simple{5}([322:409]))
% ylim([0 4])
% set(gca,'ytick',0:1:4)
% xlabel('Time')
% ylabel('Rail Selected')
% legend('Old PRs','New PRs')
% title('Pressure Rail Selections for Actuator 1 - Cap Side')
% % % 
figure(16)
subplot(2,2,1);
plot([1:321],-F1([1:321]))
xlabel('Time')
ylabel('Drive Cycle Force (N)')
title('Segment - 1 of Rail Force (No different decisions')

subplot(2,2,2)
plot([322:409],-F1([322:409]))
xlabel('Time')
ylabel('Drive Cycle Force (N)')
title('Segment - 2 of Rail Force (Different decisions)')

subplot(2,2,3)
plot([410:1296],-F1([410:1296]))
xlabel('Time')
ylabel('Drive Cycle Force (N)')
title('Segment - 3 of Rail Force (No different decisions')

subplot(2,2,4)
plot([1297:1306],-F1([1297:1306]))
xlabel('Time')
ylabel('Drive Cycle Force (N)')
title('Segment - 4 of Rail Force (Different decisions')
% 
% %% Considering HECM E_motor operating points:
% %change indices to have the required HECM plot
% figure(17)
% subplot(2,1,1);
% plot(HECM_angvel{3}, HECM_torque{3},'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]') %original
% subplot(2,1,2);
% plot(HECM_angvel{7}, HECM_torque{7},'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]') %modified

%% checking the forces at the instants where different decisions were taken:

%actuator 3
to_be_checked_forces_3 = -F3(to_be_checked_3);

%preparing histogram for actuator 3:
figure(18)
histogram(to_be_checked_forces_3,10) % splitting the histogram into 10 equal ranges
xlabel('Drive Cycle Forces')
ylabel('Number of different decisions')
title('Relation between drive cycle force and number of different decisions')

%% plotting the 3d graph that relates d.c.force and Pa and Pb:

AA = [0 0.5 1]*Pmax; % old pressure rails
BB = [0.09 0.5359 1]*Pmax; % new pressure rails

[aa_old bb_old] = ndgrid(AA,AA);
[aa_new bb_new] = ndgrid(BB,BB);

aa_old = aa_old(:);
bb_old = bb_old(:);
aa_new = aa_new(:);
bb_new = bb_new(:);

F_low_bound_old = aa_old.*ACap1;
F_up_bound_old = aa_old.*ACap1 - bb_old.*ARod1 - (MaxT1_Act*ARod1*2*pi)/(D_HECM1);

F_low_bound_new = aa_new.*ACap1;
F_up_bound_new = aa_new.*ACap1 - bb_new.*ARod1 - (MaxT1_Act*ARod1*2*pi)/(D_HECM1);

[xx_old yy_old] = meshgrid(AA', AA');
zz_up_old = reshape(F_up_bound_old,3,3);
zz_low_old = reshape(F_low_bound_old,3,3);

[xx_new yy_new] = meshgrid(BB',BB');
zz_up_new = reshape(F_up_bound_new,3,3);
zz_low_new = reshape(F_low_bound_new,3,3);

figure(19)
mesh(xx_old,yy_old,zz_low_old') %% cavitation based forces for old PR
hold on
mesh(xx_old,yy_old,zz_up_old') %% torque limit based forces for old PR
hold on
mesh(xx_new,yy_new,zz_low_new')
hold on
mesh(xx_new,yy_new,zz_up_new')
hold on
plot3(aa_old,bb_old,F_up_bound_old,'o')
hold on
plot3(aa_old,bb_old,F_low_bound_old,'o')
hold on
plot3(aa_new,bb_new,F_up_bound_new,'x')
hold on
plot3(aa_new,bb_new,F_low_bound_new,'x')
xlabel('P_A')
ylabel('P_B')
title('Plot of force limits v/s Pressure Selections')
