% Make plots which compare things from all cases and PTOS
clear, close all

%% Input Data
% Electric motor generator sizes
HHEA_main_reg_gen = 846; %kW HHEA main generator size for regular waves
HHEA_HECM_reg_gen = 1566; %kW HHEA HECM generator size for regular waves
EHA_noco_reg_gen = 3548; %kW EHA generator size for regular waves with no codesign
EHA_co_reg_gen = 2387; %kW EHA generator size for regular waves with codesign
Check_reg_gen = 862; %kW Check valve generator size for regular waves

HHEA_main_irreg_gen = 117; %kW HHEA main generator size for irregular waves
HHEA_HECM_irreg_gen = 2231; %kW HHEA HECM generator size for irregular waves
EHA_noco_irreg_gen = 4371; %kW EHA generator size for irregular waves with no codesign
EHA_co_irreg_gen = 3080; %kW EHA generator size for irregular waves with codesign
Check_irreg_gen = 1318; %kW Check valve generator size for irregular waves


%Hydraulic pump motor sizes
HHEA_main_reg_pump = 921; %cc HHEA main pump size for regular waves
HHEA_HECM_reg_pump = 5029; %cc HHEA HECM pump size for regular waves
EHA_noco_reg_pump = 5029 ; %cc EHA pump size for regular waves with no codesign
EHA_co_reg_pump = 3745; %cc EHA pump size for regular waves with codesign
Check_reg_pump = 1166; %cc Check valve pump size for regular waves

HHEA_main_irreg_pump = 262; %cc HHEA main pump size for irregular waves
HHEA_HECM_irreg_pump = 5885; %cc HHEA HECM pump size for irregular waves
EHA_noco_irreg_pump = 5885; %cc EHA pump size for irregular waves with no codesign
EHA_co_irreg_pump = 4815; %cc EHA pump size for irregular waves with codesign
Check_irreg_pump = 549.5; %cc Check valve pump size for irregular waves

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
