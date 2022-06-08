% Run optimization with battery contraint
load JCB_LongerDC.mat
% masterScale = 1.25;
masterScale = 1;
% MaxT1_Act = inf;
% MaxT2_Act = inf;
% MaxT3_Act = inf;
% MaxT4_Act = inf;

maxRPM = 5000; %RPM %used to size the HECM

MaxT1_Act = 24; 
MaxT2_Act = 78;
MaxT3_Act = 39;
MaxT4_Act = 43;

for MaxT2_Act = MaxT2_Act:78
    disp(['Now, T1 is: ',num2str(MaxT1_Act)])
    %type of drive cycle operation:
%     cycle = NinetyDeg;
%     cycle = Grading;
    cycle = Trenching;

    load 28point1ccnolock_A6VM_4Q_20191210.mat

    % Define an array of pressures
    PR = [0 1];
%     PR = [0 1/2 1]; 
    % PR = [0 1/3 2/3 1]; 
    % PR = [0 1/4 1/2 3/4 1]; 

    Make_losses_general;
    %return
    %%
    % 2 rail option
    %all_options = {'C0'}


    % % 3 rail options
    all_options = {'C00'};

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

%         figure(3);
%         subplot(221); plot(t, -F1, '-', t([1,end]),ones(2,1)*Frange1', '-', t, CPR_F1(d),'-.','linewidth',2); ylabel('Force1 N'); xlabel('Time - s');grid
%         subplot(222); plot(t, -F2, '-', t([1,end]),ones(2,1)*Frange2', '-', t, CPR_F2(d),'-.','linewidth',2); ylabel('Force2 N'); xlabel('Time - s');grid
%         subplot(223); plot(t, -F3, '-', t([1,end]),ones(2,1)*Frange3', '-', t, CPR_F3(d),'-.','linewidth',2); ylabel('Force3 N'); xlabel('Time - s');grid
%         subplot(224); plot(t, -T4, '-', t([1,end]),ones(2,1)*Trange4', '-', t, CPR_T1(d),'-.','linewidth',2); ylabel('Torque4 Nm'); xlabel('Time - s');grid

        %% Check flow/speed
%         figure(4);
%         subplot(221); plot(t, V1, '-', t, -Q_HECM_1(d_ind)/ARod1,'-.'); ylabel('Speed1 m/s'); xlabel('Time - s');grid
%         subplot(222); plot(t, V2, '-', t, -Q_HECM_2(d_ind)/ARod2,'-.'); ylabel('Speed2 m/s'); xlabel('Time - s');grid
%         subplot(223); plot(t, V3, '-', t, -Q_HECM_3(d_ind)/ARod3,'-.'); ylabel('Speed3 m/s'); xlabel('Time - s');grid
%         subplot(224); plot(t, W4, '-', t, w4(d_ind),'-.'); ylabel('Speed4 rad/s'); xlabel('Time - s');grid

        %% Check HECM E-motor operating points - (w,T)
%         figure(5)
%         subplot(211); plot(w1(d_ind), T1_Act(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]')
%         subplot(212); plot(w2(d_ind), T2_Act(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]')
%         max_power_peak_cornor = ...
%             [max(abs(w1(d_ind).*T1_Act(d_ind))),max(abs(w1(d_ind))*max(abs(T1_Act(d_ind))));
%              max(abs(w2(d_ind).*T2_Act(d_ind))),max(abs(w2(d_ind))*max(abs(T2_Act(d_ind))))]

        %% Check HECM P/M operating points - (w,P) vs HECM PM efficiency
        %Mech=T1_Act_Mapping.*w1rad_Mapping;
        %Hyd=Q1_Act_Mapping.*P1_Mapping;
        %Eff=Hyd./(Mech+eps).*(Hyd>0)+Mech./(Hyd+eps).*(Hyd<0);
        disp('Calculating HECM hydraulic efficiencies')
        HydHECM=DeltaP_HECM1(d_ind).*Q_HECM_1(d_ind);
        MechHECM=w1(d_ind).*T1_Act(d_ind);    
        indP=HydHECM>0; eff(1,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
        indM=HydHECM<0; eff(1,2)=sum(MechHECM(indM))/sum(HydHECM(indM));
        HydHECM=DeltaP_HECM2(d_ind).*Q_HECM_2(d_ind);
        MechHECM=w2(d_ind).*T2_Act(d_ind);
        indP=HydHECM>0; eff(2,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
        indM=HydHECM<0; eff(2,2)=sum(MechHECM(indM))/sum(HydHECM(indM));

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

%         figure(6); clf;
%         subplot(211); plot(w1(d_ind), DeltaP_HECM1(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]');
%         hold on; [c1,h1]=contour(w1rad_Mapping,P1_Mapping,Eff,[-0.1:0.1:1]); clabel(c1,h1)
%         axis(1.05*[min(w1(d_ind)),max(w1(d_ind)),min(DeltaP_HECM1(d_ind)),max(DeltaP_HECM1(d_ind))])
%         title(['Pumping/Motoring efficiency: ',num2str(eff(1,:)*100),'%'])
% 
%         subplot(212); plot(w2(d_ind), DeltaP_HECM2(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]')
%         hold on; [c2,h2]=contour(w1rad_Mapping,P1_Mapping,Eff,[-0.1:0.1:1]); clabel(c2,h2)
%         axis(1.05*[min(w2(d_ind)),max(w2(d_ind)),min(DeltaP_HECM2(d_ind)),max(DeltaP_HECM2(d_ind))])
%         title(['Pumping/Motoring efficiency: ',num2str(eff(2,:)*100),'%'])
    %     [TotalMPWork,NetWorkout+Totallosses]
        EnergyConservation = [Netbattery-(Totallosses+NetWorkout)]
    end
    NetRail_all
    Eff_list = Eff_all'

    [PosWork/1000 sum(ActualMPLoss)/1000 sum(HECMLosses(d_ind))*dt/1000]
    [NegWork/1000 (NetWorkout+Totallosses)/1000]

    %EnergyConservation = [Netbattery-(Totallosses+NetWorkout)]
    % sum([PosWork/1000 sum(ActualMPLoss)/1000 sum(HECMLosses(d_ind))*dt/1000]) / sum([NegWork/1000 NetWorkout+Totallosses/1000])
    
    if sum(NetRail_all(2:end)>= -100) + sum(NetRail_all(2:end) <= 100) == length(PR)*2 - 2
        disp(['MaxT1_Act = ',num2str(MaxT1_Act)])
        return
%     else 
%         disp(['MaxT2_Act = ',num2str(MaxT2_Act),' is not a solution'])
    end
end
