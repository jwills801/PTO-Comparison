% Run optimization with battery contraint

Make_losses_WEC;

Noise = rand(size(HECMLosses))*1000;
HECMLosses_noisy = HECMLosses + Noise;
%%
choices=['M','P','0']; % motoring, pumping, or nothing

% 3 rail options
all_options = {'NMM', 'NPP','NMP','NPM','C00','C0P','CP0','C0M','CM0'}

% 4 rail options
% k=0;
% for k1=1:3
%     for k2=1:3
%         for k3=1:3
%             k=k+1;
%             if (choices(k1)=='0' | choices(k2) == '0' | choices(k3) == '0'),
%                 all_options{k}=['C',choices(k1),choices(k2),choices(k3)];
%             else 
%                 all_options{k}=['N',choices(k1),choices(k2),choices(k3)];
%             end
%         end
%     end
% end

%all_options = {'C0MM'} 
tic
Eff_all=[];
NetRail_all=[];
batt_all=[];
NetRailEnergy=[];
all_options = {'Z00'}
for case_no = 1:length(all_options)
    options = all_options{case_no}
    [MPstr,Jstr,lam]=expand_options(options,PR);
    
    eval(MPstr);
    eval(['f=@(lam) -sum(',Jstr,');']);
    options = optimset('PlotFcns',@optimplotfval);
    tic
        lam = fminsearch(f,lam,options)
    toc
%% evaluate the optimal solution
    eval(['[cost,d]=',Jstr,';']) % d is the rail combo choice (#s from 1 to 9)
    
    sz = size(HECMLosses); 
    I_t = [1:sz(1)]';
    dt = mean(diff(t));
    d_ind = sub2ind(sz,I_t,d);
    battery = battery_power(d_ind);
    for kk=1:length(PR)
        NetRailEnergy(kk)=sum(QR{kk}(d_ind))*PR(kk)*dt;
    end
    NetRailEnergy;
    ActualMPLoss = NetRailEnergy.*((NetRailEnergy>0).*[0,R_P_loss'] ...
                                 - (NetRailEnergy<0).*[0,R_M_loss']);
    TotalMPWork = sum(NetRailEnergy+ActualMPLoss);
    Netbattery = sum(battery)*dt
    NetWorkin = -sum(F1.*V1)*dt;
    Totallosses = sum(HECMLosses(d_ind))*dt+sum(ActualMPLoss);
    tmp=F1.*V1; %PosWork=-sum(tmp(tmp<0))*dt;
    Efficiency = (TotalMPWork+Netbattery)/NetWorkin
    [abs(TotalMPWork+Netbattery)+Totallosses,abs(NetWorkin)]
    [sum(cost)*dt,Totallosses]

    NetRail_all(case_no,:) = NetRailEnergy
    Eff_all(case_no)=Efficiency
    batt_all(case_no) = Netbattery

    figure(2);
    subplot(211); plot(t, cumsum(battery)*dt); c_battery = cumsum(battery);
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
    figure(3);
    plot(t, -F1, '-', t([1,end]),ones(2,1)*Frange1', '-', t, CPR_F1(d),'-o','linewidth',2); 
    ylabel('Force1 N'); xlabel('Time - s');grid

    %% Check flow/speed
    figure(4);
    plot(t, V1, '-', t, -Q_HECM_1(d_ind)/ARod1,'-.'); ylabel('Speed1 m/s'); xlabel('Time - s');grid
    %% Check HECM E-motor operating points - (w,T)
    figure(5)
    plot(w1(d_ind), T1_Act(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]')
    max_power_peak_cornor = ...
        [max(abs(w1(d_ind).*T1_Act(d_ind))),max(abs(w1(d_ind))*max(abs(T1_Act(d_ind))))]

    %% Check HECM P/M operating points - (w,P) vs HECM PM efficiency
    Mech=T1_Act_Mapping.*w1rad_Mapping;
    Hyd=Q1_Act_Mapping.*P1_Mapping;
    Eff=Hyd./(Mech+eps).*(Hyd>0)+Mech./(Hyd+eps).*(Hyd<0);
    
    disp('Calculating HECM hydraulic efficiencies')
    HydHECM=DeltaP_HECM1(d_ind).*Q_HECM_1(d_ind);
    MechHECM=w1(d_ind).*T1_Act(d_ind);
    %indP=HydHECM>0; eff(1,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
    %indM=HydHECM<0; eff(1,2)=sum(MechHECM(indM))/sum(HydHECM(indM));
    indP=HydHECM>0; eff(1,1)=sum(HydHECM(indP))/sum(battery(indP));
    indM=HydHECM<0; eff(1,2)=sum(battery(indM))/sum(HydHECM(indM));
    
    figure(6); clf;
    plot(w1(d_ind), DeltaP_HECM1(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]');
    hold on; [c1,h1]=contour(w1rad_Mapping,P1_Mapping,Eff,[-0.1:0.1:1]); clabel(c1,h1)
    axis(1.05*[min(w1(d_ind)),max(w1(d_ind)),min(DeltaP_HECM1(d_ind)),max(DeltaP_HECM1(d_ind))])
    title(['Pumping/Motoring efficiency: ',num2str(eff(1,:)*100),'%'])

    pause(1)
end
NetRail_all
Eff_all'
toc
