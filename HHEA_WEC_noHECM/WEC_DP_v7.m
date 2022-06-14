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
Spectrum = 'irregular';

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
decision_ind = NaN(size(t)); % initialize decision_ind - which holds the index of the rail chosen at each time
for i = 1:length(t)
    [~,decision_ind(i)] = min( abs( f_discr - F1(i) ) ) ; % Find which force rail is being using
    cap_rail = find(PR == PRA1(decision_ind(i))) ; % Cap side rail
    rod_rail = find(PR == PRB1(decision_ind(i))) ; % Rod side rail
    
    Q(cap_rail,i) = V1(i)*ACap1 ;
    Q(rod_rail,i) = -V1(i)*ARod1 ; % velocity of piston is positive in extension AND this flow is flow leaving the accumulator
end

vol_over_time = cumsum(Q')*dt; % m^3
figure, plot(t,vol_over_time/1e3), xlabel('Time (s)'), ylabel('Cumulative Volume (L)'), legend('Tank','Low Pressure Rail','Middle Pressure Rail','High Pressure Rail','Location','Northwest')

vol = vol_over_time(end,2:end)';





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
dt = t(2); % give dt back its value after being switched by the Swtiching loss matrix code

%%
switching_loss = zeros(size(t));
for t_ind = 2:length(t)-1
    switching_loss(t_ind) = Eloss_SwitchingValve(decision_ind(t_ind-1),decision_ind(t_ind),t_ind);
end
total_switching_loss = sum(switching_loss);

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

    

    %% check results

    % checkpower
    check(1) = std( [chi_MT*D*wMP*PR(3) chi_HT*D*wMP*PR(4) chi_HL*D*wMP*(PR(4)-PR(2)) ]) < 1 ;
    
    %check time
    t_HL = vol(1)/chi_HL/D/wMP;
    t_HT = (-vol(3)-vol(1)) /chi_HT/D/wMP;
    t_MT = -vol(2)/chi_MT/D/wMP;
    check(2) = abs( t_HL+t_HT+t_MT - T) < 1e-3 ;

    %check flow contraint
    check(3) = abs(chi_HL*D*wMP*t_HL -vol(1)) < 1e-3;
    check(4) = abs(chi_MT*D*wMP*t_MT + vol(2)) < 1e-3;
    check(5) = abs(chi_HT*D*wMP*t_HT+chi_HL*D*wMP*t_HL + vol(3)) < 1e-3;
    
    if sum(check) ~= length(check)
        disp('Contraints not met - see checks')
    end

    % how efficient is the pump at each operating condition?
    Calc_Eff
    total_MP_losses = (1-total_MP_e)*chi_HT*D*wMP*PR(4)*T; % J 
else
    disp('The code here is invalid')
    return
end

gen_size = chi_HT*D*wMP*PR(4) - total_MP_losses/T - sum(switching_loss)/T;
disp(['Generator size is ', num2str(gen_size/1e3,4),' kW'])
%% Plots
% Drive cycle
figure, subplot(211), plot(t,F1/1e6), ylabel('Force (MN)')
        subplot(212), plot(t,V1), ylabel('Velocity (m/s)'), xlabel('Time (s)')

% Energy plots
WorkOut = [zeros(1,find(t==simu.rampTime)),t(2)*cumsum(gen_size*ones(1,length(t)-find(t==simu.rampTime)))];
figure, plot(t,-t(2)*cumsum(F1.*V1)/1e6,t,WorkOut/1e6), xlabel('Time (s)'), ylabel('Work (MJ)'),legend('Mechanical Energy In','Electrical Energy Out')

%% Power 
figure, plot(t,-(F1.*V1)/1e3,t,gen_size), legend('Power In','Power Out')
xlim([100 125]), ylabel('Power (kW)'), xlabel('Time (s)')