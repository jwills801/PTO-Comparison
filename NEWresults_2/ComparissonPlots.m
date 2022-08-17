% Make plots which compare things from all cases and PTOS
clear, close all

%% Input Data
% Electric motor generator sizes
HHEA_main_reg_gen = 1102; %kW HHEA main generator size for regular waves
HHEA_HECM_reg_gen = 1338; %kW HHEA HECM generator size for regular waves
EHA_noco_reg_gen = 4106; %kW EHA generator size for regular waves with no codesign
EHA_co_reg_gen = 2642; %kW EHA generator size for regular waves with codesign
Check_reg_gen = 509; %kW Check valve generator size for regular waves
Disc_HHEA_reg_gen = 1139; %kW Discrete HHEA generator size for regular waves


HHEA_main_irreg_gen = 322; %kW HHEA main generator size for irregular waves
HHEA_HECM_irreg_gen = 1301; %kW HHEA HECM generator size for irregular waves
EHA_noco_irreg_gen = 3744; %kW EHA generator size for irregular waves with no codesign
EHA_co_irreg_gen = 2685; %kW EHA generator size for irregular waves with codesign
Check_irreg_gen = 231; %kW Check valve generator size for irregular waves
Disc_HHEA_irreg_gen = 374; %kW Discrete HHEA generator size for irregular waves



%Hydraulic pump motor sizes
HHEA_main_reg_pump = 3255; %cc HHEA main pump size for regular waves
HHEA_HECM_reg_pump = 7267; %cc HHEA HECM pump size for regular waves
EHA_noco_reg_pump = 7267; %cc EHA pump size for regular waves with no codesign
EHA_co_reg_pump = 5298; %cc EHA pump size for regular waves with codesign
Check_reg_pump = 2331; %cc Check valve pump size for regular waves
Disc_HHEA_reg_pump = 1958; %kW Discrete HHEA pump size for regular waves


HHEA_main_irreg_pump = 902; %cc HHEA main pump size for irregular waves
HHEA_HECM_irreg_pump = 6972; %cc HHEA HECM pump size for irregular waves
EHA_noco_irreg_pump = 6972; %cc EHA pump size for irregular waves with no codesign
EHA_co_irreg_pump = 4916; %cc EHA pump size for irregular waves with codesign
Check_irreg_pump = 1097; %cc Check valve pump size for irregular waves
Disc_HHEA_irreg_pump = 1487; %kW Discrete HHEA pump size for irregular waves



% CWR in
HHEA_reg_CWR_in = 1.3658;
EHA_noco_reg_CWR_in = 1.3658;
EHA_co_reg_CWR_in = 1.2258;
Check_reg_CWR_in = 0.5919; 
Disc_HHEA_reg_CWR_in = 1.3809;


HHEA_irreg_CWR_in = 0.97741; 
EHA_noco_irreg_CWR_in = 0.97741;
EHA_co_irreg_CWR_in = 0.9266;
Check_irreg_CWR_in = .5408; 
Disc_HHEA_irreg_CWR_in = .9899;


% CWR out
HHEA_reg_CWR_out = 1.0255;
EHA_noco_reg_CWR_out = 0.68138;
EHA_co_reg_CWR_out = 0.83352;
Check_reg_CWR_out = 0.5070; 
Disc_HHEA_reg_CWR_out = 1.0207;


HHEA_irreg_CWR_out = 0.94194; 
EHA_noco_irreg_CWR_out = 0.46689;
EHA_co_irreg_CWR_out = 0.56671; 
Check_irreg_CWR_out = 0.4219; 
Disc_HHEA_irreg_CWR_out = .79978;
%% Bar graphs for generators
% Regular
figure
Gen_Sizes = [HHEA_main_reg_gen HHEA_HECM_reg_gen ; ... % HHEA
             EHA_noco_reg_gen  0;                  ... % EHA - no codesign
             EHA_co_reg_gen    0;                  ... % EHA - with codesign
             Check_reg_gen     0;                  ... % Check
             Disc_HHEA_reg_gen 0                   ... % Discrete HHEA
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check','Discrete HHEA'});
b = bar(labels,Gen_Sizes,'stacked');
ylabel('Generator Size (kW)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Gen_Sizes,2);
labels1 = string(Gen_Sizes(:,1)) + ' kW' + [[' + '] + string(Gen_Sizes(1,2)) + [' kW'];'';'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')
ylim([0 1.1*max(sum(Gen_Sizes,2))])

% Irregular
figure
Gen_Sizes = [HHEA_main_irreg_gen HHEA_HECM_irreg_gen ; ... % HHEA
             EHA_noco_irreg_gen  0;                  ... % EHA - no codesign
             EHA_co_irreg_gen    0;                  ... % EHA - with codesign
             Check_irreg_gen     0;                  ... % Check
             Disc_HHEA_irreg_gen 0                   ... % Discrete HHEA
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check','Discrete HHEA'});
b = bar(labels,Gen_Sizes,'stacked');
ylabel('Generator Size (kW)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Gen_Sizes,2);
labels1 = string(Gen_Sizes(:,1)) + ' kW' + [[' + '] + string(Gen_Sizes(1,2)) + [' kW'];'';'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')
ylim([0 1.1*max(sum(Gen_Sizes,2))])

%% Bar graphs for pumps
% Regular
figure
Pump_Sizes = [HHEA_main_reg_pump HHEA_HECM_reg_pump; ... % HHEA
             EHA_noco_reg_pump   0;                    ... % EHA - no codesign
             EHA_co_reg_pump     0;                    ... % EHA - with codesign
             Check_reg_pump      0;                    ... % Check
             Disc_HHEA_reg_pump  0                     ... % Discrete HHEA
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check','Discrete HHEA'});
b = bar(labels,Pump_Sizes,'stacked');
ylabel('Pump/Motor Size (cc)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Pump_Sizes,2);
labels1 = string(Pump_Sizes(:,1)) + [' cc'] + [[' + '] + string(Pump_Sizes(1,2)) + [' cc'];'';'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')
ylim([0 1.1*max(sum(Pump_Sizes,2))])

% Irregular
figure
Pump_Sizes = [HHEA_main_irreg_pump HHEA_HECM_irreg_pump; ... % HHEA
             EHA_noco_irreg_pump   0;                    ... % EHA - no codesign
             EHA_co_irreg_pump     0;                    ... % EHA - with codesign
             Check_irreg_pump      0;                    ... % Check
             Disc_HHEA_irreg_pump  0                    ... % Discrete HHEA
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check','Discrete HHEA'});
b = bar(labels,Pump_Sizes,'stacked');
ylabel('Pump/Motor Size (cc)')
xtips1 = b(1).XEndPoints;
ytips1 = sum(Pump_Sizes,2);
labels1 = string(Pump_Sizes(:,1)) + [' cc'] + [[' + '] + string(Pump_Sizes(1,2)) + [' cc'];'';'';'';''];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Main', 'HECM')
ylim([0 1.1*max(sum(Pump_Sizes,2))])

%% Bar graphs for CWR in
% Regular
figure
CWR_in = [HHEA_reg_CWR_in         HHEA_reg_CWR_out;     ... % HHEA
             EHA_noco_reg_CWR_in  EHA_noco_reg_CWR_out; ... % EHA - no codesign
             EHA_co_reg_CWR_in    EHA_co_reg_CWR_out;   ... % EHA - with codesign
             Check_reg_CWR_in     Check_reg_CWR_out;    ... % Check
             Disc_HHEA_reg_CWR_in Disc_HHEA_reg_CWR_out ... % Discrete HHEA
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check','Discrete HHEA'});
bar(labels,CWR_in);
ylabel('Capture Width Ratio')
grid
ylim([0 1.1*max(CWR_in(:))])
legend('Mechanical','Electrical','Location','Northwest')

% Irregular
figure
CWR_in = [HHEA_irreg_CWR_in         HHEA_irreg_CWR_out;      ... % HHEA
             EHA_noco_irreg_CWR_in  EHA_noco_irreg_CWR_out; ... % EHA - no codesign
             EHA_co_irreg_CWR_in    EHA_co_irreg_CWR_out;   ... % EHA - with codesign
             Check_irreg_CWR_in     Check_irreg_CWR_out;    ... % Check
             Disc_HHEA_irreg_CWR_in Disc_HHEA_irreg_CWR_out ... % Discrete HHEA
             ];
labels = categorical({'HHEA','EHA','EHA codesign','Check','Discrete HHEA'});
bar(labels,CWR_in);
ylabel('Capture Width Ratio')
grid
ylim([0 1.1*max(CWR_in(:))])
legend('Mechanical','Electrical','Location','Northwest')


% Plot all work in and out over time for each PTO
%% Regular Waves
figure, hold on, ylabel('Energy [MJ]'), xlabel('Time [s]')

% HHEA
load('HHEA_DP/Reg_HHEA_workOverTime.mat')
plot(t,Work_in_over_time/1e6,t,Work_out_over_time/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Reg_EHA_nocodesign_workOverTime.mat')
plot(t,Work_out/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Reg_EHA_codesign_workOverTime.mat')
plot(t,Work_in/1e6,t,Work_out/1e6)

% Check Valve
load('CheckValveCode/OSWEC_Hydraulic_PTO/Reg_check_WorkOverTime.mat')
plot(t,Work_mech/1e6,t,Work_elec/1e6)

% Discrete HHEA
load HHEA_noHECM/Reg_Work_over_time_Disc_HHEA.mat
plot(t,Work_in_over_time/1e6,'--')
plot(t,Work_out_over_time/1e6,'b--')


hold off, legend('Work in - no codesign', 'Work out HHEA', 'Work out EHA - no codesign', 'Work in - codesign', 'Work out EHA - codesign', 'Work in check valve', 'Work out check valve','Work in discrete HHEA','Work out discrete HHEA','Location','Northwest')


%% Irregular Waves
clear
figure, hold on, ylabel('Energy [MJ]'), xlabel('Time [s]')

% HHEA
load('HHEA_DP/Irreg_HHEA_workOverTime.mat')
plot(t,Work_in_over_time/1e6,t,Work_out_over_time/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Irreg_EHA_nocodesign_workOverTime.mat')
plot(t,Work_out/1e6)

% EHA - No codesign
load('EHA_codesign/OSWEC/Irreg_EHA_codesign_workOverTime.mat')
plot(t,Work_in/1e6,t,Work_out/1e6)

% Check Valve
load('CheckValveCode/OSWEC_Hydraulic_PTO/Irreg_check_WorkOverTime.mat')
plot(t,Work_mech/1e6,t,Work_elec/1e6)

% Discrete HHEA
load HHEA_noHECM/Irreg_Work_over_time_Disc_HHEA.mat
plot(t,Work_in_over_time/1e6,'--')
plot(t,Work_out_over_time/1e6,'b--')

hold off, legend('Work in - no codesign', 'Work out HHEA', 'Work out EHA - no codesign', 'Work in - codesign', 'Work out EHA - codesign', 'Work in check valve', 'Work out check valve','Work in discrete HHEA','Work out discrete HHEA','Location','Northwest')
