%% This will compare the different PTOs

%load('AdamsGitHub/OSWEC_Hydraulic_PTO/output_reg/wec_ptosim_out.mat')
load Jackson_CheckValve_Regular.mat
checkvalve = struct();
checkvalve.t = output.controller.time;
checkvalve.dt = mean(diff(checkvalve.t));
checkvalve.Fin = output.controller.force;
checkvalve.Vin = output.controller.velocity;
checkvalve.gentorque = output.controller.generatorTorque;
checkvalve.gen_w = output.controller.generatorAngSpeed;
checkvalve.flow = output.controller.flow;

checkvalve.Work_out = cumsum(checkvalve.gentorque.*checkvalve.gen_w)*checkvalve.dt;
checkvalve.e = -checkvalve.Work_out(end)/sum(checkvalve.Fin.*checkvalve.Vin*checkvalve.dt);
checkvalve.GenSize = max(checkvalve.gentorque.*checkvalve.gen_w); % w
checkvalve.motorsize = max(checkvalve.flow./checkvalve.gen_w)*2*pi*1e6;

figure(1)
yyaxis left
plot(checkvalve.t,-checkvalve.Fin.*checkvalve.Vin,checkvalve.t,checkvalve.gen_w.*checkvalve.gentorque)
ylim([-2.5e5/3 2.5e5])
ylabel('Power (W)')
yyaxis right
plot(checkvalve.t,checkvalve.Vin)
legend('Power in','Power out','Velocity','location','northwest')
xlim([150 200])
ylim([-1 3])
grid on
ylabel('Velocity (m/s)')
xlabel('Time (s)')
print('CheckValvePower', '-dpng', '-r600')

%%
run('../VariableDisplacement_Wec/main_Scaled.m')
MYonepass
%% 
close all
figure(2)
labels = categorical({'EHA','HHEA','Check Valve'});
b = bar(labels,[e,Efficiency,checkvalve.e]);
ylabel('Efficiency')
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylim([0 1])


figure(3)
labels = categorical({'EHA','HHEA','Check Valve'});
b = bar(labels,[EHA_Generator/1000, 0; mp_Generator/1000 HECM_Generator/1000; checkvalve.GenSize/1000 0],'stacked');
ylabel('Size of generator (kW)')
legend('Main pump','HECM','Location','northwest')
xtips1 = b(2).XEndPoints;
ytips1 = b(2).YEndPoints;
labels1 = string(round(ytips1)) + [' kW'];
labels1(2) = string(round(mp_Generator/1000)) + [' kW +'] +  string(round(HECM_Generator/1000)) + [' kW'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title(['Electric Generator Sizing ( Total Mean Power = ',num2str(-Total_Work_OverTime(end)/t(end)/1000),' kW)'])
ylim([0 1000])

figure(4)
labels = categorical({'EHA','HHEA','Check Valve'});
b = bar(labels,[EHA_Pump_size, 0; MP_size HECM_size; checkvalve.motorsize 0],'stacked');
ylabel('Size of motor (cc)')
legend('Main pump/motor','HECM pump/motor','Location','northwest')
xtips1 = b(2).XEndPoints;
ytips1 = b(2).YEndPoints;
labels1 = string(round(ytips1)) +  [' cc'];
labels1(2) = string(round(MP_size)) + [' cc +'] +  string(round(HECM_size)) + [' cc'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Hydraulic Pump/Motor Sizing')
ylim([0 1.25*(MP_size + HECM_size)])

figure(5)
plot(t, -F1, '-', t([1,end]),ones(2,1)*Frange1', '-', t, CPR_F1(d),'-o','linewidth',2); 
ylabel('Force (N)'); xlabel('Time (s)');grid
title('HHEA rail selections')

figure(6)
P1 = plot(t,cumsum(-P_in)*dt,t,cumsum(-P_out)*dt,t,-Total_Work_OverTime,checkvalve.t,-cumsum(checkvalve.Fin.*checkvalve.Vin)*checkvalve.dt,checkvalve.t,checkvalve.Work_out);
legend('Work in controlled','Electrical Energy out EHA','Electrical Energy out HHEA','Work in Check Valve','Electrical Energy Out Check Valve')
ylabel('Work (J)')
xlabel('Time (s)')
xlim([0 200])
title('Energy Absorbtion')

