clear, close all
%% TO DO:
% Rerun wec sim with finer sampling time
%%
disp('____________________________________________________________________________')
%cycle='regular';
cycle='irregular';

rail_no = 3;
%rail_no = 4;

eval(['load ../PI_',cycle,'_',num2str(rail_no),'Rail','.mat'])

Emap = 'constantEff';

V1 = myoutput.signals.values(11,:)';
F1 = -myoutput.signals.values(17,:)';
t = myoutput.time;

DPdt = t(2); % Seconds between switches
t_c = 0:DPdt:t(end);
spacing = round(DPdt/t(2)); % re-sampling rate - we can make decisions every "spacing" milliseconds

if round(t_c(end),4) ~= round(t(end),4)
    error('Constrained time does not fill up the whole space')
end
maxRPM = 2000; %RPM

% Define an array of pressures
Pmax = 25e6; % I've been using the same Pmax for reg and irreg
if rail_no == 3
    PR = [0 2/3 1]*Pmax;
elseif rail_no == 4
    PR = [0, 1/6, 2/3, 1]*Pmax;
end

Make_losses_DP;

GetSwitchingLossMatrix_v3;
%Eloss_A1 = 0*Eloss_A1;
%Eloss_B1 = 0*Eloss_B1;

Sum_to_constrain_switches_v1;

NetWorkout = -sum(F1.*V1)*t(2);
OutputPowerSum = sum(abs(F1.*V1))*t(2);

all_options = gen_options(length(PR), ['0','M','P']);
all_options = {'CMM'}; % Best for REGULAR, 3 rails
%all_options = {''}; % Best for IRREGULAR

for case_no = 1:length(all_options)
    option = all_options{case_no}
    [Jstr,lam]=expand_options(option,PR);
    %HECMLosses_c = HECMLosses_c + 0*rand(size(HECMLosses_c));
    eval(['f=@(lam) ',Jstr,';']);


    params = struct(); params.Eloss_A1 = Eloss_A1; params.Eloss_B1 = Eloss_B1;
    params.t_c = t_c; params.f = f; params.PR = PR; params.PRA1 = PRA1; params.PRB1 = PRB1;
    
    %options = optimset('PlotFcns',@optimplotfval,'MaxIter',100);
    options = optimset('MaxIter',100);
    tic
    [lam,cost1,exitflag] = fminsearch(@(lam) -fun(lam,params),lam,options);
    Convergence = exitflag==1;
    disp(['fminsearch took ' num2str(toc) ' seconds'])

    %% evaluate the optimal solution
    GetDecisions

    battery = battery_power(d_ind);
    for kk=1:length(PR)
        NetRailEnergy(kk)=sum(QR{kk}(d_ind))*PR(kk)*t(2);
    end
    ActualMPLoss = NetRailEnergy.*((NetRailEnergy>0).*[0,R_P_loss'] ...
        - (NetRailEnergy<0).*[0,R_M_loss']);
    TotalMPWork = sum(NetRailEnergy+ActualMPLoss);
    HECM_hyd_loss = sum(PLoss_Hydraulic_HECM(d_ind))*t(2);
    HECM_elect_loss = sum(PLoss_Electric_HECM(d_ind))*t(2);
    HECM_total_loss = sum(HECMLosses(d_ind))*t(2);
    MP_total_loss = sum(ActualMPLoss);

    switching_loss = zeros(size(t_c));
    d_vector_c = [starting_ind, Decision_vector_c]; 
    if Tf > t_c(2)
        switch_spacing = Tf/t_c(2);
    else
        switch_spacing = 1;
    end
    t_ind_c_prev = 2;
    for t_ind_c = 2:switch_spacing:length(t_c)
        if 1%d_vector_c(t_ind_c_prev) ~= d_vector_c(t_ind_c)
            ind_A1_old = find(PRA1(d_vector_c(t_ind_c_prev)) == PR); ind_B1_old = find(PRB1(d_vector_c(t_ind_c_prev)) == PR);
            ind_A1_new = find(PRA1(d_vector_c(t_ind_c)) == PR); ind_B1_new = find(PRB1(d_vector_c(t_ind_c)) == PR);
            switching_loss(t_ind_c-1) = Eloss_A1(ind_A1_old, ind_A1_new, t_ind_c ) + Eloss_B1(ind_B1_old, ind_B1_new, t_ind_c );
        end
        t_ind_c_prev = t_ind_c;
    end
    TotalSwitchingLoss = sum(switching_loss); % No dt because each entry is already a energy unit, not a power
    Netbattery = sum(battery)*t(2);
    Totallosses = HECM_total_loss+MP_total_loss+TotalSwitchingLoss;
    tmp=F1.*V1; PosWork=-sum(tmp(tmp<0))*t(2);
    NegWork = PosWork - NetWorkout;
    Efficiency = (TotalMPWork+TotalSwitchingLoss + Netbattery)/NetWorkout;
    corner = max(abs(w1(d_ind))*max(abs(T1_Act(d_ind))));
    disp('Optimized cost vs Total losses')
    disp([cost*t_c(2),Totallosses]) % if these don't match, it could be that the option was not consistent
    % cost matches total losses if you account for the lam*NetrailEnergy and the NetRailEnergy*R_P/M_Loss parts when the option is '0'
    disp('Energy input vs energy output and losses')
    disp([TotalMPWork+Netbattery+NegWork+TotalSwitchingLoss,PosWork+MP_total_loss+HECM_total_loss+TotalSwitchingLoss]) % note the switching loss is both an input and a loss

    corner_all(case_no)=corner;
    NetRail_all(case_no,:) = NetRailEnergy;
    Eff_all(case_no)=Efficiency;
end
if length(Eff_all) > 1
    disp(Eff_all)
    [a,b] = max(Eff_all);   
    disp(['Best option is number ', num2str(b), ' - ', all_options{b} ' - which is ' num2str(round(a*100,1)) '% efficient'])
else
%%
figure;
%subplot(211); plot(t, cumsum(battery)*t(2)); ylabel('J'), c_battery = cumsum(battery)*t(2);
%title(['Battery end state = ',num2str(c_battery(end)),'J'])
if  length(PR)==5
    %subplot(212); 
    plot(t, cumsum([QR{1}(d_ind)',QR{2}(d_ind)',QR{3}(d_ind)',QR{4}(d_ind)',QR{5}(d_ind)'])*t(2));
    legend('QR{1}','QR{2}','QR{3}','QR{4}','QR{5}');
elseif length(PR)==4
    %subplot(212); 
    plot(t, cumsum([QR{1}(d_ind)',QR{2}(d_ind)',QR{3}(d_ind)',QR{4}(d_ind)'])*t(2));
    legend('Tank',[num2str(round(PR(2)/1e6,2)) ' MPa rail'],[num2str(round(PR(3)/1e6,2)) ' MPa rail'],[num2str(round(PR(4)/1e6,2)) ' MPa rail'],'Location','NorthWest');
elseif length(PR)==3
    %subplot(212); 
    plot(t, cumsum([QR{1}(d_ind)',QR{2}(d_ind)',QR{3}(d_ind)'])*t(2));
    legend('QR{1}','QR{2}','QR{3}');
end
ylabel('Cumulative Flow [m^3]');
xlabel('Time [s]');
%end
%save SC_5rails_results.mat NetRail_all Eff_all PR ScaleMaxT optimal
%toc

%% Plot drive cycle force and velocity
figure, yyaxis left, plot(t,-F1/1e6), ylabel('Force [MN]'), xlabel('Time [s]'), ylim([-10 10])
yyaxis right, plot(t,V1), ylabel('Velocity [m/s]'), xlim([100 125])

%% Check force/torque
%figure, plot(t, -F1, '-', t([1,end]),ones(2,1)*Frange1', '-', t_c, CPR_F1(d_vector_c),'x','linewidth',2); ylabel('Force1 N'); xlabel('Time - s');grid
figure(111), hold on, plot(t_c, CPR_F1(d_vector_c)/1e6,'x','MarkerSize',12,'linewidth',2), hold off
xlim([100, 125])

%% Switching losses
figure, plot(t,switching_loss*1e-3), ylabel('Switching Loss [kJ]'), xlabel('Time [s]'), grid on, xlim([100 125])

%% Check flow/speed
%figure,plot(t, V1, '-', t, -Q_HECM_1/ARod1,'-.'); ylabel('Speed1 m/s'); xlabel('Time - s');grid

%% plot torque
figure, plot(t, T1_Act(d_ind), '-'); ylabel('Torque [Nm]'); xlabel('Time  [s]'); grid, hold on
plot([0 t(end)], [MaxT1_Act MaxT1_Act],'--k'), plot([0 t(end)], [-MaxT1_Act -MaxT1_Act],'--k'), hold off, ylim(1.25*[-MaxT1_Act MaxT1_Act]), xlim([0 t(end)])
xlim([100 125])

%% Plot Flow
%figure, plot(t, Q_HECM_1, '-'); ylabel('Flow 1 (m^3/s)'); xlabel('Time - s');grid

%% Plot Pressure Drop
%figure(7), plot(t, DeltaP_HECM1(d_ind)/1e6, '-'); ylabel('Delta P 1 (MPa)'); xlabel('Time - s');grid

%% Plot Pressure at Rod Side - outlet of HECM
%figure, plot(t, P1_Cyl(d_ind)/1e6, '-'); ylabel('Rod Pressure 1 (MPa)'); xlabel('Time - s');grid,title('Pressure at Outlet of HECM 1')%, xlim([100 130])

%% Energy Distribution
figure
labels = categorical({'HECM Energy','Main Pump Energy'});
b = bar(labels,[-Netbattery/1e6, -TotalMPWork/1e6-TotalSwitchingLoss/1e6]);
ylabel('Work (MJ)')
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData)+ [' MJ'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([min([0 -Netbattery/1e6]) -1.1*TotalMPWork/1e6])

%% Loss Distribution
Totallosses = HECM_total_loss+MP_total_loss+TotalSwitchingLoss;
figure
labels = categorical({'HECM Losses','Main Pump Losses','Throttling Losses'});
b = bar(labels,[HECM_total_loss,MP_total_loss, TotalSwitchingLoss]/1e6);
ylabel('Energy [MJ]')
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(round(b(1).YData,1))+ [' MJ'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
%ylim([0 13])

%%
t_start = 50;
disp('Losses:'), disp(['   HECM:      ',num2str(HECM_total_loss,3)]), disp(['   Main Pump: ',num2str(MP_total_loss,3)]), disp(['   Switching: ',num2str(TotalSwitchingLoss,3)])
disp(['Waves: ' cycle])

disp(['HECM hydraulic unit is ' num2str(D_HECM1*1e6),' cc'])
size_MP;
disp(['MP hydraulic unit is ' num2str(D_MP*1e6),' cc'])

[max_batt_elec,max_batt_elec_ind] = max(abs(battery));
[max_batt_mech,max_batt_mech_ind] = max(abs(T1_Act(d_ind).*w1(d_ind)));
if max_batt_elec > max_batt_mech
    disp(['Size of HECM generator: ' num2str(max_batt_elec/1000) ' kW, set at time ' num2str(t(max_batt_elec_ind)) ' seconds.'])
else
    disp(['Size of HECM generator: ' num2str(max_batt_mech/1000) ' kW, set at time ' num2str(t(max_batt_mech_ind)) ' seconds.'])
end
disp(['Size of MP generator: ' num2str(-(TotalMPWork + TotalSwitchingLoss)/(t(end)-t_start)/1000/.9) ' kW'])
disp(['Overall Efficiency: ' num2str(Efficiency) ''])

Work_in_over_time = cumsum(F1.*V1)*t(2);
Work_out_over_time = [zeros(1,round(t_start/t(2))+1), (-(TotalMPWork + TotalSwitchingLoss)/(t(end)-t_start))*t(2)*(0:1:length(t)-(round(t_start/t(2))+2))]-cumsum(battery)*t(2);
figure, plot(t,Work_in_over_time/1e6,t,Work_out_over_time/1e6), ylabel('Work [MJ]'), xlabel('Time [s]')

figure, plot(t,F1.*V1/1e6,t([1 round(t_start/t(2)) round(t_start/t(2))+1 end]), [0 0 -(TotalMPWork + TotalSwitchingLoss)/(t(end)-t_start) -(TotalMPWork + TotalSwitchingLoss)/(t(end)-t_start)]/1e6)
legend('Power in','Power out'), ylabel('Power [MW]'), xlabel('Time [s]'), xlim([100 125])

% Find Capture Wdith Ratio
B = 18; % m (This is the width of the oswec)
ave_power_in = (Work_in_over_time(end) - Work_in_over_time(round(t_start/t(2))+1))/(t(end)-t_start); % Ave power for last 150 s
CW_in = ave_power_in/waves.Pw;
disp(['CWR in: ' num2str(CW_in/B) ])
ave_power_out = (Work_out_over_time(end) - Work_out_over_time(round(t_start/t(2))+1))/(t(end)-t_start); % Ave power for last 150 s
CW_out = ave_power_out/waves.Pw;
disp(['CWR out: ' num2str(CW_out/B) ])

end


function cost = fun(lam,params)

% unpack params
t_c = params.t_c; f = params.f; PR = params.PR;
PRA1 = params.PRA1; PRB1 = params.PRB1;
Eloss_A1 = params.Eloss_A1; Eloss_B1 = params.Eloss_B1;

CostMatrix = f(lam);

% The actual function
J = NaN(length(t_c),length(PRA1));
J(end,:) = zeros(size(PRA1));
for t_ind = (length(t_c)-1):-1:1
    maxSwitchingLoss = max(Eloss_A1(:, :, t_ind+1 ),[],'All') + max(Eloss_B1(:, :, t_ind+1 ),[],'All');    
    [min_w_out_switches, min_ind_w_out_switches] = min(CostMatrix(t_ind+1,:) + J(t_ind+1,:));
    inds_2_check = find(   (CostMatrix(t_ind+1,:) + J(t_ind+1,:))  < maxSwitchingLoss/t_c(2) + min_w_out_switches );

    SwitchingLoss = (PRA1==PR)*Eloss_A1(:, (PRA1(min_ind_w_out_switches)==PR), t_ind+1 ) + (PRB1==PR)*Eloss_B1(:, (PRB1(min_ind_w_out_switches)==PR), t_ind+1 );

    J(t_ind,:) = min_w_out_switches + SwitchingLoss/t_c(2);
    
    for i = inds_2_check

        % Construct Switching loss vector for this specific time index and
        % valve configuration
        % i denotes the previous valve configuration, we are deciding where
        % to go next by taking the min 
        SwitchingLoss = Eloss_A1((PR==PRA1(i)), :, t_ind+1 )*(PRA1==PR)' + Eloss_B1((PR==PRB1(i)), :, t_ind+1 )*(PRB1==PR)';

        J(t_ind,i) = min( CostMatrix(t_ind+1,:) + J(t_ind+1,:) + SwitchingLoss/t_c(2) );
        % ELoss is divided by the time step because it is an energy, while
        % the other terms are power
    end
end

cost = min(J(1,:));


end % fun
