% Run optimization with battery contraint
%close all
MyMake_losses_WEC;

Noise = rand(size(HECMLosses))*1000;
HECMLosses_noisy = HECMLosses + Noise;
%%
choices=['M','P','0']; % motoring, pumping, or nothing
solutionsFound = [];
eff = [];

% 3 rail options
%all_options = {'NMM', 'NPP','NMP','NPM','C00','C0P','CP0','C0M','CM0'};

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
% %all_options = {'C0MM'} 

Eff_all=[];
NetRail_all=[];
batt_all=[];
NetRailEnergy=[];
all_options = {'Z00M'}; % Z00M for 4 PR (ScaleTmax = 2.3), Z0M for 3PR (ScaleTmax = 1.1)

for case_no = 1:length(all_options)
    options = all_options{case_no}
    [MPstr,Jstr,lam]=MYexpand_options(options,PR);
    
    eval(MPstr);
    eval(['f=@(lam) -sum(',Jstr,');']);
    options = optimset('PlotFcns',@optimplotfval);   
    lam = fminsearch(f,lam,options);
    
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
    ActualMPLoss = NetRailEnergy.*((NetRailEnergy>0).*[0,R_P_loss'] ...
                                 - (NetRailEnergy<0).*[0,R_M_loss']);
    TotalMPWork = sum(NetRailEnergy+ActualMPLoss);
    Netbattery = sum(battery)*dt;
    NetWorkin = -sum(F1.*V1)*dt;
    %%
    Totallosses = sum(HECMLosses(d_ind))*dt+sum(ActualMPLoss);
    tmp=F1.*V1; %PosWork=-sum(tmp(tmp<0))*dt;
    Efficiency = (TotalMPWork+Netbattery)/NetWorkin
    [(TotalMPWork+Netbattery)-Totallosses,(NetWorkin)];
    [sum(cost)*dt,Totallosses];


    % Energy from MP is spread out over the whole time period
    MPstaringt = 50; % seconds (must be a multiple of dt, otherwise the MP_work line of code needs to be changed
    MP_work = [zeros(find(t==MPstaringt),1); cumsum(TotalMPWork/(t(end)-MPstaringt)*ones(length(t)-find(t==MPstaringt),1))*dt];
    HECM_work = cumsum(battery)*dt;
    Total_Work_OverTime = MP_work + HECM_work;
    MP_power = [zeros(find(t==MPstaringt),1); (TotalMPWork/(t(end)-MPstaringt)*ones(length(t)-find(t==MPstaringt),1))];
    
    
    NetRail_all(case_no,:) = NetRailEnergy;
    Eff_all(case_no)=Efficiency;
    batt_all(case_no) = Netbattery;

%     figure(2);
%     subplot(211); plot(t, -cumsum(battery)*dt); c_battery = cumsum(battery);
%     title(['Battery end state = ',num2str(c_battery(end)),'J'])
%     if  length(PR)==5
%         subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind),QR{5}(d_ind)])*dt); 
%         legend('QR{1}','QR{2}','QR{3}','QR{4}','QR{5}');
%     elseif length(PR)==4
%         subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind)])*dt); 
%         legend('QR{1}','QR{2}','QR{3}','QR{4}');
%     elseif length(PR)==3
%         subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind)])*dt); 
%         legend('QR{1}','QR{2}','QR{3}');
%     elseif length(PR)==2
%         subplot(212); plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind)])*dt); 
%         legend('QR{1}','QR{2}');
%     end
%     ylabel('Cumulative Rail Flow');
%     xlabel('Time - s');
%end
%save SC_5rails_results.mat NetRail_all Eff_all PR ScaleMaxT optimal
%toc
    %% Check force/torque
    figure;
    plot(t, -F1, '-',t([1,end]),ones(2,1)*Frange1', '-','linewidth',2)
    hold on
    scatter( t, CPR_F1(d),'o','MarkerEdgeColor',[0.6350 0.0780 0.1840]); 
    hold off
    ylabel('Force (N)'); xlabel('Time (s)');grid
    xlim([150 200])
    if length(PR) == 4
        print('Figures/4RailsDecided', '-dpng', '-r600')
    elseif length(PR) == 3
        print('Figures/3RailsDecided', '-dpng', '-r600')
    else
    end
    % title('4 Pressure Rails')

    %% Show torque of HECM during that same time interval
        figure(33);
    plot(t,T1_Act(d_ind), '-','linewidth',2)
    ylabel('Tourque (Nm)'); xlabel('Time (s)');grid
    xlim([150 200])
    print('Figures/3Rails_HECM_TOURQUE', '-dpng', '-r600')
    % title('4 Pressure Rails')
    

    %% Rail Flow
    figure
    grid on
    if  length(PR)==5
        plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind),QR{5}(d_ind)])*dt); 
        legend('QR{1}','QR{2}','QR{3}','QR{4}','QR{5}');
    elseif length(PR)==4
        plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind)])*dt); 
        legend([num2str(PR(1)*1e-6) ' MPa Pressure Rail'],[num2str(PR(2)*1e-6) ' MPa Pressure Rail'],[num2str(PR(3)*1e-6) ' MPa Pressure Rail'],[num2str(PR(4)*1e-6) ' MPa Pressure Rail'],'location','northwest');
        print('Figures/4RailFlow', '-dpng', '-r600')
    elseif length(PR)==3
        plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind)])*dt); 
        legend([num2str(PR(1)*1e-6) ' MPa Pressure Rail'],[num2str(PR(2)*1e-6,3) ' MPa Pressure Rail'],[num2str(PR(3)*1e-6,3) ' MPa Pressure Rail'],'location','northwest');
        print('Figures/3RailFlow', '-dpng', '-r600')
    elseif length(PR)==2
        plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind)])*dt); 
        legend('QR{1}','QR{2}');
    end
    ylabel('Cumulative Rail Flow ($m^3$)','interpreter','latex');
    xlabel('Time  ($s$)','interpreter','latex');

    
    
    %% Check flow/speed
%     figure(4);
%     plot(t, V1, '-', t, -Q_HECM_1(d_ind)/ARod1,'-.'); ylabel('Speed1 m/s'); xlabel('Time - s');grid
%     %% Check HECM E-motor operating points - (w,T)
%     figure(5)
%     plot(w1(d_ind), T1_Act(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Torque [Nm]')
%     max_power_peak_cornor = ...
%         [max(abs(w1(d_ind).*T1_Act(d_ind))),max(abs(w1(d_ind))*max(abs(T1_Act(d_ind))))]
% 
%     %% Check HECM P/M operating points - (w,P) vs HECM PM efficiency
%     Mech=T1_Act_Mapping.*w1rad_Mapping;
%     Hyd=Q1_Act_Mapping.*P1_Mapping;
%     Eff=Hyd./(Mech+eps).*(Hyd>0)+Mech./(Hyd+eps).*(Hyd<0);
%     
%     disp('Calculating HECM hydraulic efficiencies')
%     HydHECM=DeltaP_HECM1(d_ind).*Q_HECM_1(d_ind);
%     MechHECM=w1(d_ind).*T1_Act(d_ind);
%     %indP=HydHECM>0; eff(1,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
%     %indM=HydHECM<0; eff(1,2)=sum(MechHECM(indM))/sum(HydHECM(indM));
%     indP=HydHECM>0; eff(1,1)=sum(HydHECM(indP))/sum(battery(indP));
%     indM=HydHECM<0; eff(1,2)=sum(battery(indM))/sum(HydHECM(indM));
%     
%     figure(6); clf;
%     plot(w1(d_ind), DeltaP_HECM1(d_ind),'.');grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]');
%     hold on; [c1,h1]=contour(w1rad_Mapping,P1_Mapping,Eff,[-0.1:0.1:1]); clabel(c1,h1)
%     axis(1.05*[min(w1(d_ind)),max(w1(d_ind)),min(DeltaP_HECM1(d_ind)),max(DeltaP_HECM1(d_ind))])
%     title(['Pumping/Motoring efficiency: ',num2str(eff(1,:)*100),'%'])

%% Check to to see if signs are consistent
% Four rail case
checker = 0;
for k=2:length(PR)
    switch all_options{case_no}(k)
        case '0'
            if abs(NetRailEnergy(k)) <= 6e5
                checker = checker + 1;
            end
        case 'M'
            if NetRailEnergy(k) <= 0
                checker = checker + 1;
            end
        case 'P'
            if NetRailEnergy(k) >= 0
                checker = checker + 1;
            end
    end
end
if checker == length(PR)-1
    all_options{case_no}
    solutionsFound = [solutionsFound;all_options{case_no}]
end
case_no
eff = [eff, Efficiency];
end
%%

% calculateds the hydraulic MP size is options are NOM
% if all_options{1} == 'Z0M'
%     V3 = sum(QR{3}(d_ind)*dt);
%     MP_size = -V3/wMP/t(end)*2*pi*1e6; % cc
%     % 75 cc - compared to 1170 cc for the EHA case
%     
%     % That is roughly 15 X smaller
%     % But the HECM motor is 749 cc
% end 

% figure(55)
% subplot(211);
% plot(t,Q_HECM_1(d_ind)./w1(d_ind)*2*pi*1e6);
% xlim([0 50])
% subplot(212)
% plot(t, -F1, '-', t([1,end]),ones(2,1)*Frange1', '-', t, CPR_F1(d),'-o','linewidth',2); 
% ylabel('Force1 N'); xlabel('Time - s');grid
% xlim([0 50])


% To size the Main hydraulic pump/motor, we will operate the motor so that its displacement changes at each 
% Pressure rail so that the main generator will be able to operate at
% constant power
ActualRailEnergy = sign(NetRailEnergy).*(abs(NetRailEnergy)- ActualMPLoss);
if length(PR) == 4
    Vol_rail2 = sum(QR{2}(d_ind)*dt);% Volume needed to be pumped / motored at rail 2
    Vol_rail3 = sum(QR{3}(d_ind)*dt);
    Vol_rail4 = sum(QR{4}(d_ind)*dt);
    
    %|P2|=|P3|=|P4|
    fracdisp2 = sign(ActualRailEnergy(2))*1;
    fracdisp3 = sign(ActualRailEnergy(3))*PR(2)/PR(3);
    fracdisp4 = sign(ActualRailEnergy(4))*PR(2)/PR(4);
    
    % t2 = t3 = t4
    MP_size = (ActualRailEnergy(2)/PR(2)/wMP/fracdisp2 + ActualRailEnergy(3)/PR(3)/wMP/fracdisp3 + ActualRailEnergy(4)/PR(4)/wMP/fracdisp4 ) /t(end) *2*pi*1e6;
    
    t2 = ActualRailEnergy(2)/PR(2)/wMP/fracdisp2/MP_size*(2*pi*1e6);
    t3 = ActualRailEnergy(3)/PR(3)/wMP/fracdisp3/MP_size*(2*pi*1e6);
    t4 = ActualRailEnergy(4)/PR(4)/wMP/fracdisp4/MP_size*(2*pi*1e6);
    
   %Check
   % See how much power is required at each time interval
   % How much power does it take to pump the second rail?
   Power2 = PR(2)*wMP*fracdisp2*MP_size/(2*pi*1e6); % Pressure times Flow rate
   Power3 = PR(3)*wMP*fracdisp3*MP_size/(2*pi*1e6);
   Power4 = PR(4)*wMP*fracdisp4*MP_size/(2*pi*1e6);
   [Power2 Power3 Power4];
   
   % Can these fracdisplacements and this MP size do the required amount of
   % flow to each rail?
   [Vol_rail2 wMP*fracdisp2*MP_size/(2*pi*1e6)*t2] % Rail 2
   [Vol_rail3 wMP*fracdisp3*MP_size/(2*pi*1e6)*t3] % Rail 3
   [Vol_rail4 wMP*fracdisp4*MP_size/(2*pi*1e6)*t4] % Rail 4
    
   [Power2*t2 ActualRailEnergy(2)]
   [Power3*t3 ActualRailEnergy(3)]
   [Power4*t4 ActualRailEnergy(4)]
   
   [t2+t3+t4 200]
   mp_Generator = Power2
   
   
   
   %% This is for the C00M case
    MP_size = abs(ActualRailEnergy(4))/PR(4)/wMP /t(end) *2*pi*1e6
    mp_Generator = PR(4)*wMP*MP_size/(2*pi*1e6)
else
    % length(PR) = 3
%     MP_size = abs(ActualRailEnergy(3))/PR(3)/wMP /t(end) *2*pi*1e6
%     mp_Generator = PR(3)*wMP*MP_size/(2*pi*1e6)
    Vol_rail2 = sum(QR{2}(d_ind)*dt);% Volume needed to be pumped / motored at rail 2
    Vol_rail3 = sum(QR{3}(d_ind)*dt);
    
    %|P2|=|P3|=|P4|
    fracdisp2 = sign(ActualRailEnergy(2))*1;
    fracdisp3 = sign(ActualRailEnergy(3))*PR(2)/PR(3);
    
    % t2 = t3 = t4
    MP_size = (ActualRailEnergy(2)/PR(2)/wMP/fracdisp2 + ActualRailEnergy(3)/PR(3)/wMP/fracdisp3 ) /t(end) *2*pi*1e6;
    
    t2 = ActualRailEnergy(2)/PR(2)/wMP/fracdisp2/MP_size*(2*pi*1e6);
    t3 = ActualRailEnergy(3)/PR(3)/wMP/fracdisp3/MP_size*(2*pi*1e6);
    
   %Check
   % See how much power is required at each time interval
   % How much power does it take to pump the second rail?
   Power2 = PR(2)*wMP*fracdisp2*MP_size/(2*pi*1e6) % Pressure times Flow rate
   Power3 = PR(3)*wMP*fracdisp3*MP_size/(2*pi*1e6)
   [Power2 Power3 Power4]
   
   % Can these fracdisplacements and this MP size do the required amount of
   % flow to each rail?
   [Vol_rail2 wMP*fracdisp2*MP_size/(2*pi*1e6)*t2] % Rail 2
   [Vol_rail3 wMP*fracdisp3*MP_size/(2*pi*1e6)*t3] % Rail 3
    
   [Power2*t2 ActualRailEnergy(2)]
   [Power3*t3 ActualRailEnergy(3)]
   
   [t2+t3+t4 200]
   mp_Generator = Power2

end




HECM_size = ScaleHECM1*107 % cc
%% Size the generators
%HECM_Generator = max(abs(battery)); % dont do this bc in the EHA case, I
%sized it by looking at the energy before the generator not after it.
HECM_Generator = max(abs(T1_Act(d_ind).*w1(d_ind))) % W


% figure(10)
% plot(t,T1_Act(d_ind).*w1(d_ind))

% The HECM generator is 750 kW
% The MP generator is 115 kW
% The EHA generator is 929 kW


% figure(11)
% plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind)])*dt); 
% legend('QR{1}','QR{2}','QR{3}','QR{4}');
% ylabel('Cumulative Rail Flow (m^3)');
% xlabel('Time (s)');
% title('Net Flows')

% figure(12)
% plot(t, (PR(3)*QR{3}(d_ind))); 
% legend('QR{3}');
% ylabel('Rail Power (W) at each time');
% xlabel('Time (s)');


%% Account for switching losses
% d(end+1) = d(end);
% d_ind(end+1) = d_ind(end);
if length(PR) ==3
SwitchingLosses = 0;
for i = 1:length(d)-1 
    %SwitchingLosses(i) = ((PRA1(d(i)) == PR(1) & PRA1(d(i+1)) == PR(2))|(PRA1(d(i)) == PR(2) & PRA1(d(i+1)) == PR(3))).*(2192684.7248*Q_HECM_1(d_ind(i)).^2+97636.3858*Q_HECM_1(d_ind(i))+5588.8859) + (PRA1(d(i)) == PR(1) & PRA1(d(i+1)) == PR(3)).*(2384217.9498*Q_HECM_1(d_ind(i)).^2+219488.4184*Q_HECM_1(d_ind(i))+19486.4537) + ((PRA1(d(i)) == PR(2) & PRA1(d(i+1)) == PR(1))|(PRA1(d(i)) == PR(3) & PRA1(d(i+1)) == PR(2))).*(2192684.7248*Q_HECM_1(d_ind(i)).^2+-97636.3858*Q_HECM_1(d_ind(i))+5588.8859) +(PRA1(d(i)) == PR(3) & PRA1(d(i+1)) == PR(1)).*(2384217.9498*Q_HECM_1(d_ind(i)).^2+-219488.4184*Q_HECM_1(d_ind(i))+19486.4537) + ...
       % ((PRB1(d(i)) == PR(1) & PRB1(d(i+1)) == PR(2))|(PRB1(d(i)) == PR(2) & PRB1(d(i+1)) == PR(3))).*(2192684.7248*Q_HECM_1(d_ind(i)).^2+97636.3858*Q_HECM_1(d_ind(i))+5588.8859) + (PRB1(d(i)) == PR(1) & PRB1(d(i+1)) == PR(3)).*(2384217.9498*Q_HECM_1(d_ind(i)).^2+219488.4184*Q_HECM_1(d_ind(i))+19486.4537) + ((PRB1(d(i)) == PR(2) & PRB1(d(i+1)) == PR(1))|(PRB1(d(i)) == PR(3) & PRB1(d(i+1)) == PR(2))).*(2192684.7248*Q_HECM_1(d_ind(i)).^2+-97636.3858*Q_HECM_1(d_ind(i))+5588.8859) +(PRB1(d(i)) == PR(3) & PRB1(d(i+1)) == PR(1)).*(2384217.9498*Q_HECM_1(d_ind(i)).^2+-219488.4184*Q_HECM_1(d_ind(i))+19486.4537);
       PRA = PRA1(d(i));
       PRB = PRA1(d(i+1));
       Q = Q_HECM_1(d_ind(i));
       P_Low = PR(1); P_Mid = PR(2); P_High = PR(3);
       SwitchingLosses(i) = ((PRA == P_Low & PRB == P_Mid)|(PRA == P_Mid & PRB == P_High)).*(1080223.9146*Q.^2+96847.0548*Q+10672.7225) + (PRA == P_Low & PRB == P_High).*(1259010.3959*Q.^2+242155.5914*Q+36677.9714) + ((PRA == P_Mid & PRB == P_Low)|(PRA == P_High & PRB == P_Mid)).*(1080223.9146*Q.^2+-96847.0548*Q+10672.7225) +(PRA == P_High & PRB == P_Low).*(1259010.3959*Q.^2+-242155.5914*Q+36677.9714);

end 
% d = d(1:end-1);
% d_ind = d_ind(1:end-1);
SwitchingLosses(end+1) = 0;
% figure(10)
% plot(Q_HECM_1(d_ind),SwitchingLosses(1:end))
% figure(11)
% plot(t,SwitchingLosses(1:end))
TotalSwitchingLosses = sum(SwitchingLosses)*dt;

% Account for switching losses in net totals. Note that energy in and out
% are negative and the losses are posative
Totallosses = Totallosses + TotalSwitchingLosses;
Efficiency = (TotalMPWork+Netbattery+TotalSwitchingLosses)/NetWorkin;
[(TotalMPWork+Netbattery+TotalSwitchingLosses)-Totallosses,(NetWorkin)];
[sum(cost)*dt,Totallosses];
Total_Work_OverTime = Total_Work_OverTime + SwitchingLosses';
abs(TotalSwitchingLosses/NetWorkin);

else
    SwitchingLosses = 0;

    P_Low = PR(1);
    P_Mid1 = PR(2);
    P_Mid2 = PR(3);
    P_High = PR(4);
    for i = 1:length(d)-1
        PRa = PRA1(d(i));
        PRb = PRA1(d(i+1));
        Q = Q_HECM_1(d_ind(i));
        SwitchingLosses(i) = (PRa == P_Low & PRb == P_Mid1).*(1149143.7229*Q.^2+141714.3799*Q+17734.2079)+ (PRa == P_Low & PRb == P_Mid2).*(885299.0598*Q.^2+23110.8247*Q+1599.4304)+(PRa == P_Low & PRb == P_High).*(1259010.3959*Q.^2+242155.5914*Q+36677.9714) + (PRa == P_Mid1 & PRb == P_Low).*(1149143.7229*Q.^2+-141714.3799*Q+17734.2079) + (PRa == P_Mid2 & PRb == P_Low).*(885299.0598*Q.^2+-23110.8247*Q+1599.4304) + (PRa == P_Mid1 & PRb == P_High).*(995905.9904*Q.^2+56697.9819*Q+5268.361) + (PRa == P_Mid2 & PRb == P_High).*(1207801.0438*Q.^2+190348.6152*Q+26411.3157) + (PRa == P_High & PRb == P_Low).*(1259010.3959*Q.^2+-242155.5914*Q+36677.9714) + (PRa == P_High & PRb == P_Mid1).*(995905.9904*Q.^2+-56697.9819*Q+5268.361) + (PRa == P_High & PRb == P_Mid2).*(1207801.0438*Q.^2+-190348.6152*Q+26411.3157);
    end
    % d = d(1:end-1);
    % d_ind = d_ind(1:end-1);
    SwitchingLosses(end+1) = 0;
    % figure(10)
    % plot(Q_HECM_1(d_ind),SwitchingLosses(1:end))
    % figure(11)
    % plot(t,SwitchingLosses(1:end))
    TotalSwitchingLosses = sum(SwitchingLosses)*dt;
    
    % Account for switching losses in net totals. Note that energy in and out
    % are negative and the losses are posative
    Totallosses = Totallosses + TotalSwitchingLosses;
    Efficiency = (TotalMPWork+Netbattery+TotalSwitchingLosses)/NetWorkin
    [(TotalMPWork+Netbattery+TotalSwitchingLosses)-Totallosses,(NetWorkin)];
    [sum(cost)*dt,Totallosses];
    Total_Work_OverTime = Total_Work_OverTime + SwitchingLosses';
    abs(TotalSwitchingLosses/NetWorkin);
    
    
end

eff'
solutionsFound
totalGen = mp_Generator+ HECM_Generator
