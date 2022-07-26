% Make plots which compare things from all cases and PTOS
clear, close all

%% Input Data
% Electric motor generator sizes
HHEA_main_reg_gen = 1220; %kW HHEA main generator size for regular waves
HHEA_HECM_reg_gen = 1454; %kW HHEA HECM generator size for regular waves
EHA_noco_reg_gen = 4098; %kW EHA generator size for regular waves with no codesign
EHA_co_reg_gen = 2886; %kW EHA generator size for regular waves with codesign

Check_reg_gen = 509; %kW Check valve generator size for regular waves


HHEA_main_irreg_gen = 329; %kW HHEA main generator size for irregular waves
HHEA_HECM_irreg_gen = 1367; %kW HHEA HECM generator size for irregular waves
EHA_noco_irreg_gen = 3739; %kW EHA generator size for irregular waves with no codesign

EHA_co_irreg_gen = 3080; %kW EHA generator size for irregular waves with codesign
Check_irreg_gen = 231; %kW Check valve generator size for irregular waves


%Hydraulic pump motor sizes
HHEA_main_reg_pump = 3438; %cc HHEA main pump size for regular waves
HHEA_HECM_reg_pump = 7461; %cc HHEA HECM pump size for regular waves
EHA_noco_reg_pump = 7461; %cc EHA pump size for regular waves with no codesign
EHA_co_reg_pump = 5885; %cc EHA pump size for regular waves with codesign

Check_reg_pump = 2331; %cc Check valve pump size for regular waves


HHEA_main_irreg_pump = 968; %cc HHEA main pump size for irregular waves
HHEA_HECM_irreg_pump = 6842; %cc HHEA HECM pump size for irregular waves
EHA_noco_irreg_pump = 6842; %cc EHA pump size for irregular waves with no codesign

EHA_co_irreg_pump = 4815; %cc EHA pump size for irregular waves with codesign
Check_irreg_pump = 1097; %cc Check valve pump size for irregular waves


% CWR out
HHEA_reg_CWR = 0.99764;
EHA_noco_reg_CWR = 0.68467;
EHA_co_reg_CWR = 0.82974;
Check_reg_CWR = ; 


HHEA_irreg_CWR = 0.93257; 
EHA_noco_irreg_CWR = 0.47242;
EHA_co_irreg_CWR = ; 
Check_irreg_CWR = ; 

%% Bar graphs for generators
% Regular
figure
Gen_Sizes = [HHEA_main_reg_gen HHEA_HECM_reg_gen ; ... % HHEA
             EHA_noco_reg_gen  0;                  ... % EHA - no codesign
             EHA_co_reg_gen    0;                  ... % EHA - with codesign
             Check_reg_gen     0                   ... % Check
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check'});
b = bar(labels,Gen_Sizes,'stacked');
ylabel('Generator Size (kW)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Gen_Sizes,2);
labels1 = string(Gen_Sizes(:,1)) + [' kW'] + [[' + '] + string(Gen_Sizes(1,2)) + [' kW'];'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')

% Irregular
figure
Gen_Sizes = [HHEA_main_irreg_gen HHEA_HECM_irreg_gen ; ... % HHEA
             EHA_noco_irreg_gen  0;                  ... % EHA - no codesign
             EHA_co_irreg_gen    0;                  ... % EHA - with codesign
             Check_irreg_gen     0                   ... % Check
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check'});
b = bar(labels,Gen_Sizes,'stacked');
ylabel('Generator Size (kW)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Gen_Sizes,2);
labels1 = string(Gen_Sizes(:,1)) + [' kW'] + [[' + '] + string(Gen_Sizes(1,2)) + [' kW'];'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')

%% Bar graphs for pumps
% Regular
figure
Pump_Sizes = [HHEA_main_reg_pump HHEA_HECM_reg_pump; ... % HHEA
             EHA_noco_reg_pump   0;                    ... % EHA - no codesign
             EHA_co_reg_pump     0;                    ... % EHA - with codesign
             Check_reg_pump      0                     ... % Check
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check'});
b = bar(labels,Pump_Sizes,'stacked');
ylabel('Pump/Motor Size (cc)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Pump_Sizes,2);
labels1 = string(Pump_Sizes(:,1)) + [' cc'] + [[' + '] + string(Pump_Sizes(1,2)) + [' cc'];'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')

% Irregular
figure
Pump_Sizes = [HHEA_main_irreg_pump HHEA_HECM_irreg_pump; ... % HHEA
             EHA_noco_irreg_pump   0;                    ... % EHA - no codesign
             EHA_co_irreg_pump     0;                    ... % EHA - with codesign
             Check_irreg_pump      0                     ... % Check
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check'});
b = bar(labels,Pump_Sizes,'stacked');
ylabel('Pump/Motor Size (cc)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Pump_Sizes,2);
labels1 = string(Pump_Sizes(:,1)) + [' cc'] + [[' + '] + string(Pump_Sizes(1,2)) + [' cc'];'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')



% Plot all work in and out over time for each PTO
%% Regular Waves
figure, hold on, ylabel('Energy (MJ)'), xlabel('Time (s)')

% HHEA
load('HHEA_DP/Reg_HHEA_workOverTime.mat')
plot(t,-WorkIn/1e6,t,-WorkOut/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Reg_EHA_nocodesign_workOverTime.mat')
plot(t,Work_out/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Reg_EHA_codesign_workOverTime.mat')
plot(t,Work_in/1e6,t,Work_out/1e6)

% Check Valve
load('CheckValveCode/OSWEC_Hydraulic_PTO/Reg_check_WorkOverTime.mat')
plot(t,Work_in/1e6,t,Work_out/1e6)

hold off, legend('Work in - no codesign', 'Work out HHEA', 'Work out EHA - no codesign', 'Work in - codesign', 'Work out EHA - codesign', 'Work in check valve', 'Work out check valve')


%% Irregular Waves
clear
figure, hold on, ylabel('Energy (MJ)'), xlabel('Time (s)')

% HHEA
load('HHEA_DP/Irreg_HHEA_workOverTime.mat')
plot(t,Work_in/1e6,t,-WorkOut/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Irreg_EHA_nocodesign_workOverTime.mat')
plot(t,Work_out/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Irreg_EHA_codesign_workOverTime.mat')
plot(t,Work_in/1e6,t,Work_out/1e6)

% Check Valve
load('CheckValveCode/OSWEC_Hydraulic_PTO/Irreg_check_WorkOverTime.mat')
plot(t,Work_in/1e6,t,Work_out/1e6)

hold off, legend('Work in - no codesign', 'Work out HHEA', 'Work out EHA - no codesign', 'Work in - codesign', 'Work out EHA - codesign', 'Work in check valve', 'Work out check valve')
