%% This code will cmake plots which compare the perforance of the check valve PTO, the EHA PTO, HHEA with 3 PR and HHEA with 4 PR
% First compare everything with 3 PR

% SET UP MYonepass.m with 3 PR stuff
% SET UP MYmake_Losses_WEC.m with 3 PR stuff


run('Compare.m')


% Save 3 PR stuff
Efficiency3 = Efficiency;
mp_Generator3 = mp_Generator;
HECM_Generator3 = HECM_Generator;
MP_size3 = MP_size;
HECM_size3 = HECM_size;
Total_Work_OverTime3 = Total_Work_OverTime;
battery3 = battery;
ActualMPLoss3 = sum(ActualMPLoss);
HECMLosses3 = sum(HECMLosses(d_ind))*dt;
TotalSwitchingLosses3 = TotalSwitchingLosses;
Totallosses3 = Totallosses;



%%
% SET UP MYonepass.m with 4 PR stuff
% SET UP MYmake_Losses_WEC.m with 4 PR stuff

run('Compare.m')


% Save 4 PR stuff
Efficiency4 = Efficiency;
mp_Generator4 = mp_Generator;
HECM_Generator4 = HECM_Generator;
MP_size4 = MP_size;
HECM_size4 = HECM_size;
Total_Work_OverTime4 = Total_Work_OverTime;
battery4 = battery;
ActualMPLoss4 = sum(ActualMPLoss);
HECMLosses4 = sum(HECMLosses(d_ind))*dt;
TotalSwitchingLosses4 = TotalSwitchingLosses;
Totallosses4 = Totallosses;


%%
%set(fig.Children, ...
    % 'FontName',     'Times', ...
    % 'FontSize',     9);
close all
figure(2)
labels = categorical({'EHA','HHEA 3 Rails','HHEA 4 Rails','Check Valve'});
b = bar(labels,[e,Efficiency3,Efficiency4,checkvalve.e]);
ylabel('Efficiency')
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0 1])
print('Results_Efficiency', '-dpng', '-r600')
%%

figure(3)
labels = categorical({'EHA','HHEA 3 Rails','HHEA 4 Rails','Check Valve'});
b = bar(labels,[EHA_Generator/1000, 0; mp_Generator3/1000 HECM_Generator3/1000;mp_Generator4/1000 HECM_Generator4/1000; checkvalve.GenSize/1000 0],'stacked');
ylabel('Size of generator (kW)')
legend('Main pump','HECM','Location','northeast')
xtips1 = b(2).XEndPoints;
ytips1 = b(2).YEndPoints;
labels1 = string(round(ytips1)) + [' kW'];
labels1(2) = string(round(mp_Generator3/1000)) + [' kW +'] +  string(round(HECM_Generator3/1000)) + [' kW'];
labels1(3) = string(round(mp_Generator4/1000)) + [' kW +'] +  string(round(HECM_Generator4/1000)) + [' kW'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
% title(['Electric Generator Sizing ( Total Mean Power = ',num2str(-Total_Work_OverTime(end)/t(end)/1000),' kW)'])
ylim([0 1000])
print('Results_Generator', '-dpng', '-r600')
%%
figure(4)
labels = categorical({'EHA','HHEA 3 Rails','HHEA 4 Rails','Check Valve'});
b = bar(labels,[EHA_Pump_size, 0; MP_size3 HECM_size3; MP_size4 HECM_size4; checkvalve.motorsize 0],'stacked');
ylabel('Size of motor (cc)')
legend('Main pump/motor','HECM pump/motor','Location','northwest')
xtips1 = b(2).XEndPoints;
ytips1 = b(2).YEndPoints;
labels1 = string(round(ytips1)) +  [' cc'];
labels1(2) = string(round(MP_size3)) + [' cc +'] +  string(round(HECM_size3)) + [' cc'];
labels1(3) = string(round(MP_size4)) + [' cc +'] +  string(round(HECM_size4)) + [' cc'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
%title('Hydraulic Pump/Motor Sizing')
ylim([0 1.25*(MP_size4 + HECM_size4)])
print('Results_Pump', '-dpng', '-r600')
%%

figure(5)
P1 = plot(t,cumsum(-P_in)*dt/1000,t,cumsum(-P_out)*dt/1000,t,-Total_Work_OverTime3/1000,t,-Total_Work_OverTime4/1000,checkvalve.t,-cumsum(checkvalve.Fin.*checkvalve.Vin)*checkvalve.dt/1000,checkvalve.t,checkvalve.Work_out/1000);
legend('Work in with PI controller','Electrical Energy out EHA','Electrical Energy out HHEA 3 Rails','Electrical Energy out HHEA 4 Rails','Work in uncontrolled','Electrical Energy Out Check Valve','Location','northwest')
ylabel('Work (kJ)')
xlabel('Time (s)')
xlim([0 200])
print('Results_PowerComparison', '-dpng', '-r600')
%title('Energy Absorbtion')

%%
figure(6)
yyaxis left
plot(t,-P_in,t,ones(size(t))*mp_Generator3,t,battery3,'k')
xlabel('Time (s)')
ylabel('Power (W)')
ylim([-3e5 12e5])
yyaxis right
plot(t,v,'Color','#D95319')
legend('Power In','Power Out Main Pump','Power Out HECM','Velocity','location','northwest')
xlim([150 200])
ylim([-3 12])
ylabel('Velocity (m/s)')
grid on


figure(7)
yyaxis left
plot(t,-P_in,t,ones(size(t))*mp_Generator4,t,battery4,'k')
xlabel('Time (s)')
ylabel('Power (W)')
ylim([-2e5 12e5])
yyaxis right
plot(t,v,'Color','#D95319')
legend('Power In','Power Out Main Pump','Power Out HECM','Velocity','location','northwest')
xlim([150 200])
ylim([-2 12])
ylabel('Velocity (m/s)')
grid on
%print('HHEAPower', '-dpng', '-r600')
%%
figure(2)
labels = categorical({'HHEA 3 Rails','HHEA 4 Rails'});
b = bar(labels,[ActualMPLoss3,HECMLosses3, TotalSwitchingLosses3; ActualMPLoss4, HECMLosses4, TotalSwitchingLosses4],'stacked');
ylabel('Energy Loss (J)')
legend('Main Pump Losses','HECM Losses','Switching Losses','Location','northwest')

b(3).FaceColor = [0.4660 0.6740 0.1880]
ylim([0 3.5e7])