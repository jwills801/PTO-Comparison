clear, close all
%The HECM is on the Rod Side, but area ratio is one
% Results heavily depend on the sizing of the switching valve.

%PR_spread = [0 .5 1];
PR_spread = [0 1/6 2/3 1];
choices=['M','P','0']; % motoring, pumping, or nothing
k = 0;
for k1=1:3
    for k2=1:3
        for k3=1:3
            k=k+1;
            if (choices(k1)=='0' | choices(k2) == '0' | choices(k3) == '0'),
                all_options{k}=['C',choices(k1),choices(k2),choices(k3)];
            else 
                all_options{k}=['N',choices(k1),choices(k2),choices(k3)];
            end
        end
    end
end


Spectrum = 'regular';
%Spectrum = 'irregular';

IncludeSwitchingLosses = 1; % 1 for yes, include switching losses, 0 for no, dont inlcude


%PM = 'Minal';
PM = '107cc';


switch Spectrum
    case 'regular'
        load('PI_regular_discretized.mat')
    case 'irregular'
        load('PI_irregular_discretized.mat')
    otherwise
        disp('Please provide a vaid spectrum type (irregular or regular)')
        pause()
end
V1 = squeeze(myoutput.signals.values(11,:,:));
F1 = squeeze(myoutput.signals.values(17,:,:));
t = myoutput.time;
dt = mean(diff(t));
t_step = dt;
ACap1 = .2382;
ARod1 = ACap1;

%all_options = {'C00M'} % maybe best case for regular discrete
all_options = all_options(21); % maybe best case for regular discrete
ScaleMaxT1 = 2.25;
ScaleMaxPow = .05;
Pmax = 35e6; % MPa

% Find Capture Wdith Ratio
B = 18; % m (This is the width of the oswec)
t_start = 50;
Energy_first_chunk = -sum(F1(1:find(t==t_start)).*V1(1:find(t==t_start)))*dt;
Total_energy_in = -sum(F1.*V1)*dt;
ave_power_out = (Total_energy_in-Energy_first_chunk)/(t(end)-t_start); % Ave power for last 100 s
CW = ave_power_out/waves.Pw;
CWR = CW/B




%% Load in Pump Maps
% MP map
load("AP_107cc_Sept16.mat")
MP_Pmap = P1_Mapping;
MP_Qmap = Q1_Act_Mapping;
MP_Tmap = T1_Act_Mapping;
MP_wmap = w1rad_Mapping;


%% Make Pressure Rails
PR = Pmax*PR_spread; 

[PRCap1, PRRod1] = ndgrid(PR,PR);
PRA1 = PRCap1(:); PRB1 = PRRod1(:); % A is cap

figure, plot(t,F1), hold on, plot(t([1,end]),repmat(PRA1*ACap1-PRB1*ARod1,1,2)), hold off, ylabel('Force (N)'), xlabel('Time (s)'), xlim([100 150])

%% Get net flow through each rail
f_discr = PRA1*ACap1-PRB1*ARod1; % Force rail options

Q = zeros(length(PR),length(t)); % Initialize flow matrix
for i = 1:length(t)
    [~,ind] = min( abs( f_discr - F1(i) ) ) ; % Find which force rail is being using
    cap_rail = find(PR == PRA1(ind)) ; % Cap side rail
    rod_rail = find(PR == PRB1(ind)) ; % Rod side rail
    
    Q(cap_rail,i) = V1(i)*ACap1 ;
    Q(rod_rail,i) = -V1(i)*ARod1 ; % velocity of piston is positive in extension AND this flow is flow leaving the accumulator
end

vol_over_time = cumsum(Q')*dt; % m^3
figure, plot(t,vol_over_time/1e3), xlabel('Time (s)'), ylabel('Cumulative Volume (L)'), legend('Tank','Low Pressure Rail','Middle Pressure Rail','High Pressure Rail','Location','Northwest')

vol = vol_over_time(end,2:end)'

%% Main pump losses
if vol(1)>0 && vol(2)<0 && vol(3) < 0 % This code only works for a specific case of main pump flows
    % Let the flow from the second rail 
    wMP = 3000*2*pi/60; % Angular Velocity [rad/s]
    
    % Find fractional displacements
    % subscripts chi_ij denote that flow is taken from rail i and delivered to rail j
    chi_MT = 1; % Full displacement at lowest pressure difference
    chi_HT = PR(3)/PR(4);
    chi_HL = PR(3)/(PR(4)-PR(2));

    % find required displacement of pump - from time and flow requirements
    T = 150 ; % How long is the main pump running?
    D = ( vol(1)/chi_HL + ( -vol(3) - vol(1)  )/chi_HT - vol(2)/chi_MT) / wMP / T ; % m^3/rad
    disp(['Pump size is ', num2str(D*2*pi*1e6,4),' cc'])

    % Find power
    [chi_MT*D*wMP*PR(3) chi_HT*D*wMP*PR(4) chi_HL*D*wMP*(PR(4)-PR(2))] / 1e3
    disp(['Generator size is ', num2str(chi_HT*D*wMP*PR(4)/1e3,4),' kW'])

    % check results
    %check time
    t_HL = vol(1)/chi_HL/D/wMP
    t_HT = (-vol(3)-vol(1)) /chi_HT/D/wMP
    t_MT = -vol(2)/chi_MT/D/wMP
    [ t_HL+t_HT+t_MT T]

    %check flow contraint
    [chi_HL*D*wMP*t_HL vol(1)]
    [chi_MT*D*wMP*t_MT -vol(2)]
    [chi_HT*D*wMP*t_HT+chi_HL*D*wMP*t_HL -vol(3)]


    % how efficient is the pump at each operating condition?
    Calc_Eff
else
    disp('The code here is invalid')
end

return
%%

% Main Pump Losses
wMP = 3000*2*pi/60; % Angular Velocity [rad/s]
T_PR = interp2(MP_Pmap,MP_wmap,MP_Tmap,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);
Q_PR = interp2(MP_Pmap,MP_wmap,MP_Qmap,[PR(2:end)';-PR(2:end)'],ones(2*(length(PR)-1),1)*wMP);

MP_PumpEff = (Q_PR(1:(length(PR)-1)).*PR(2:end)')./(T_PR(1:(length(PR)-1))*wMP);
MP_MotorEff = T_PR(length(PR):end)*wMP./(-PR(2:end)'.*Q_PR(length(PR):end));

R_P_loss = (1./MP_PumpEff-1);  % loss per power
R_M_loss = (1-MP_MotorEff);

HECM_Q1 = repmat(-V1'*ARod1,length(PR)^2,1); 
for k=1:length(PR)
    QR{k} = zeros(length(PRA1),length(t));
    ind = find(PRA1==PR(k)); QR{k}(ind,:) = QR{k}(ind,:)+ACap1*ones(length(ind),1)*V1';
    ind = find(PRB1==PR(k)); QR{k}(ind,:) = QR{k}(ind,:)+HECM_Q1(ind,:);
end


%% Find Switching Loss Matrix
tic
Vol_inHoses = 2e-1; % m^3
Vol_A1 = Vol_inHoses; Vol_B1 = Vol_inHoses;

k_ = max(abs(V1))*ACap1/sqrt(5e5); % Q/sqrt(delP) Q --> max rated Q for the valve
zeta = 0.7; %damping coefficient
wn = 50*2*pi; % 50 Hz, 50*2*pi rad/s
beta = 1.8e9; %bulk modulus %pure oil - 1.8, typical oil mixture - approx 1.5

Tau = 1/zeta/wn; % Time constant
Tf = round(10*Tau,2);
dt = 1e-5;
tspan = 0:dt:Tf;
gs = tf(wn^2,[1,2*wn*zeta,wn^2]);

maxdelay = round(4.1*Tau/dt);
delayvals = linspace(1,maxdelay,5);

if IncludeSwitchingLosses == 1
GetSwitchingLossMatrix_v3; % for 3 rails 5 delayvals, and dt = 1e-5, it takes .5 hour
                              % for 4 rails 5 delayvals, and dt = 1e-5, it takes .8 hour
else
    Eloss_SwitchingValve = zeros(length(PR)^4,length(PR)^4,length(t));
end

toc



%% Analyze results
Net_work_in = -sum(F1.*V1)*t(2);
tmp=F1.*V1; PosWork=sum(tmp(tmp>0))*(t(2)-t(1));
NegWork = sum(tmp(tmp<0))*(t(2)-t(1));


d_battery_power = battery_power(d_ind);
d_battery_work = cumsum(d_battery_power)*(t(2)-t(1));
for kk=1:length(PR)
    NetRailEnergy(kk)=sum(QR{kk}(d_ind))*PR(kk)*(t(2)-t(1));
end
ActualMPLoss = NetRailEnergy.*((NetRailEnergy>0).*[0,R_P_loss'] - (NetRailEnergy<0).*[0,R_M_loss']);
TotalMPWork = sum(NetRailEnergy+ActualMPLoss);
Netbattery = sum(d_battery_power)*(t(2)-t(1));

d_Switchingloss = NaN(size(t)); d_Switchingloss(1) = 0; d_Switchingloss(end) = 0;
for t_ind = 2:length(t)-1
    d_Switchingloss(t_ind) = Eloss_SwitchingValve(Decision_vector(t_ind-1),Decision_vector(t_ind),t_ind+1);
end
Total_input_Energy = TotalMPWork + d_battery_work(end) + sum(d_Switchingloss);
Efficiency(case_no) = -Total_input_Energy/Net_work_in;
Actual_cost(case_no) = sum(ActualMPLoss)/t_step + sum(HECM_PowLoss(d_ind)) +sum(d_Switchingloss)/t_step;

%figure, plot(t,d_Switchingloss)

%% Make plots
pie22=[PosWork, sum(ActualMPLoss), sum(HECM_PowLoss(d_ind))*(t(2)-t(1)), IncludeSwitchingLosses*sum(d_Switchingloss)];
Totallosses = sum(HECM_PowLoss(d_ind))*(t(2)-t(1)) +sum(ActualMPLoss) + IncludeSwitchingLosses*sum(d_Switchingloss); %Joules
[PosWork+NegWork+Totallosses, Total_input_Energy]
% Battery Power
%figure, plot(t,d_battery_power), ylabel('Power (W)'), xlabel('Time (s)'), title('Battery Power'), grid on

% Battery work
figure, plot(t,d_battery_work), ylabel('Work (J)'), xlabel('Time (s)'), title(['Net Battery Work = ',num2str(d_battery_work(end)),' J']), grid on

% Net flow and main pump size
D_MP = 0;
figure, hold on
for i = 1:length(PR)
    eval(['plot(t,cumsum(QR{',num2str(i),'}(d_ind))*(t(2)-t(1)))'])
    eval(['D_MP = D_MP + sum(QR{',num2str(i),'}(d_ind))*(t(2)-t(1)) * PR(',num2str(i),')/PR(2) ']);
end, hold off, legend, grid on, ylabel('Cumulative Flow(m^3)'), xlabel('Time (s)')
D_MP = D_MP/150/2000*60*1e6 %cc Main pump size

%% Generator sizes
HECM_gen = max(abs(d_battery_power))/1e3 % kW
MP_gen = abs(TotalMPWork)/t(end)/1e3 % kW



% Force rail section
figure, plot(t,F1), hold on, plot(t([1,end]),repmat(PRA1*ACap1-PRB1*ARod1,1,2)),scatter(t,PRA1([starting_ind Decision_vector])*ACap1-PRB1([starting_ind Decision_vector])*ARod1,'x'), hold off, ylabel('Force (N)'),xlim([100 200]), xlabel('Time (s)')

% Drive cycle
figure, subplot(211), plot(t,F1/1e6), ylabel('Force (MN)')
        subplot(212), plot(t,V1), ylabel('Velocity (m/s)'), xlabel('Time (s)')

% Torque
figure, plot(t,HECM_T1(d_ind)), ylabel('Torque (Nm)'), xlabel('Time (s)'), grid on, xlim([100 150]), xlabel('Time (s)')
%hold on, plot([t(1) t(end)],[MaxT1_Act, MaxT1_Act],'--'), hold off

% Angular Speed
figure, plot(t,HECM_w1(d_ind)), ylabel('Angular Speed (rad/s)'), xlabel('Time (s)'), title('HECM Angular Speed'), grid on, xlim([100 200]), xlabel('Time (s)')

% Energy plots
MPstart = 50; %Seconds
WorkOut = d_battery_work + [zeros(1,find(t==MPstart)),t(2)*cumsum(TotalMPWork/(t(end)-MPstart)*ones(1,length(t)-find(t==MPstart)))] + cumsum(d_Switchingloss');
figure, plot(t,-t(2)*cumsum(F1.*V1)/1e6,t,-WorkOut/1e6), xlabel('Time (s)'), ylabel('Work (MJ)'),legend('Mechanical Energy In','Electrical Energy Out')

%% Power 
figure, plot(t,-(F1.*V1)/1e3,t,[zeros(1,find(t==MPstart)),-TotalMPWork/(t(end)-MPstart)*ones(1,length(t)-find(t==MPstart))]/1e3,'--k',t,(d_battery_power+d_Switchingloss'/t(2))/1e3,'--r'), legend('Power In','Main Pump Power Out','HECM Power Out')
xlim([100 125]), ylabel('Power (kW)'), xlabel('Time (s)')
%% Energy Distribution
figure
labels = categorical({'HECM Energy','Main Pump Energy'});
vol = bar(labels,[-d_battery_work(end)/1e6, -TotalMPWork/1e6]);
ylabel('Work (MJ)')
xtips1 = vol(1).XEndPoints;
ytips1 = vol(1).YEndPoints;
labels1 = string(vol(1).YData)+ [' MJ'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0 -1.1*TotalMPWork/1e6])

%% Loss Distribution
figure
Totallosses = sum(HECM_PowLoss(d_ind))*(t(2)-t(1)) +sum(ActualMPLoss) + IncludeSwitchingLosses*sum(d_Switchingloss); %Joules
labels = categorical({'HECM Losses','Main Pump Losses','Throttling Losses'});
vol = bar(labels,[sum(HECM_PowLoss(d_ind))*(t(2)-t(1))/1e6, sum(ActualMPLoss)/1e6, sum(d_Switchingloss)/1e6]);
ylabel('Energy (MJ)')
xtips1 = vol(1).XEndPoints;
ytips1 = vol(1).YEndPoints;
labels1 = string(vol(1).YData)+ [' MJ'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')



%% Maximum acceleration of HECM
% d_HECM_w1 = HECM_w1(d_ind); d_HECM_w2 = HECM_w2(d_ind);
% [a_w1,b_w1] = max(diff(d_HECM_w1)/t(2)); [a_w2,b_w2] = max(diff(d_HECM_w2)/t(2));
% [a_w,b_w] = max([a_w1,a_w2]);
% if b_w == 1
%     disp(['Maximum angular acceleration is on the lift acuator and occurs at time t = ',num2str(b_w1*t(2)),' seconds'])
% elseif b_w == 2
%     disp(['Maximum angular acceleration is on the tilt acuator and occurs at time t = ',num2str(b_w2*t(2)),' seconds'])
% else
%     disp('Something went wrong....'), pause()
% end

%% Calculating HECM hydraulic efficiencies
HydHECM=HECM_DeltaP1(d_ind).*(HECM_Q1(d_ind));
MechHECM=HECM_w1(d_ind).*HECM_T1(d_ind);
indP=HydHECM>0; eff(1,1)=sum(HydHECM(indP))/sum(MechHECM(indP));
indM=HydHECM<0; eff(1,2)=sum(MechHECM(indM))/sum(HydHECM(indM));


switch PM
    case '107cc'
        figure, plot(HECM_w1(d_ind), HECM_DeltaP1(d_ind),'.') ;grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]');
        hold on; [c1,h1]=contour(w1rad_Mapping,P1_Mapping,ElectricEff*Eff,[-0.1:0.1:1],'LineWidth',3); clabel(c1,h1,'FontSize',15)
        axis(1.05*[min(HECM_w1(d_ind)),max(HECM_w1(d_ind)),min(HECM_DeltaP1(d_ind)),max(HECM_DeltaP1(d_ind))])
        title(['Pumping/Motoring efficiency: ',num2str(eff(1,:)*100),'%'])
    case 'Minal'
        figure, plot(HECM_w1(d_ind), HECM_DeltaP1(d_ind),'.') ;grid; xlabel('Speed [rad/s]'); ylabel('Pressure [Pa]');
        hold on; [c1,h1]=contour(repmat(w_x,length(P_y),1),repmat(P_y,1,length(w_x)),Total_Efficiency,[-0.1:0.1:1],'LineWidth',3); clabel(c1,h1,'FontSize',15)
        axis(1.05*[min(HECM_w1(d_ind)),max(HECM_w1(d_ind)),(max(min(HECM_DeltaP1(d_ind)),min(P_y))),min(max(HECM_DeltaP1(d_ind)),max(P_y))])
        title(['Pumping/Motoring efficiency: ',num2str(eff(1,:)*100),'%'])
    otherwise
        disp('Please definie HECM P/M map to use')
        pause;
end



%% Save.mat file
save('LossMatrices.mat','HECM_PowLoss','QR','R_M_loss','R_P_loss','PR','Eloss_SwitchingValve')



%% Objective Function
function cost = fun(lam,params)

% unpack params
t = params.t; PRA1 = params.PRA1; f = params.f; IncludeSwitchingLosses = params.IncludeSwitchingLosses; Eloss_SwitchingValve = params.Eloss_SwitchingValve;

CostMatrix = f(lam);

% The actual function
J = NaN(length(PRA1),length(t));
J(:,length(t)) = zeros(size(PRA1));
for t_ind = (length(t)-1):-1:1
    for i = 1:length(PRA1)
        if IncludeSwitchingLosses == 0
            J(i,t_ind) = min( CostMatrix(:,t_ind+1) + J(:,t_ind+1) );
        elseif IncludeSwitchingLosses == 1
            J(i,t_ind) = min( CostMatrix(:,t_ind+1) + J(:,t_ind+1) + transpose(Eloss_SwitchingValve(i,:,t_ind+1))/(t(2)-t(1)) );
        else
            disp('Please provide a boolean value for "IncludeSwithingLosses"')
            pause()
            
        end
    end
end

cost = min(J(:,1));


end