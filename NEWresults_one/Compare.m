%% This will compare the different PTOs

load Jackson_check_four.mat
checkvalve = struct();
checkvalve.t = output.controller.time;
checkvalve.dt = output.controller.time(2) - output.controller.time(1);
checkvalve.F = -output.controller.force;
checkvalve.V = output.controller.velocity;
checkvalve.T = output.controller.torque_generator;
checkvalve.w = output.controller.angvel_generator;
checkvalve.e = sum(checkvalve.T.*checkvalve.w)/sum(checkvalve.F.*checkvalve.V);
checkvalve.Work_in_overtime = cumsum(checkvalve.F.*checkvalve.V)*checkvalve.dt;
checkvalve.Work_out_overtime = cumsum(checkvalve.T.*checkvalve.w)*checkvalve.dt;
checkvalve.GenSize = max(abs(checkvalve.T.*checkvalve.w));
checkvalve.motorsize= max(abs(checkvalve.V))*ptosim.pistonNCF.topA*60/3000/2/pi*1e6;


close all
figure
yyaxis left
plot(checkvalve.t,checkvalve.F.*checkvalve.V/1000,checkvalve.t,checkvalve.w.*checkvalve.T/1000)
ylim([-20 1000])
ylabel('Power (lW)')
yyaxis right
plot(checkvalve.t,checkvalve.V)
legend('Power in','Power out','Velocity','location','northwest')
xlim([150 180])
ylim([-1 7])
grid on
ylabel('Velocity (m/s)')
xlabel('Time (s)')
print('Figures/CheckValvePower', '-dpng', '-r600')

figure()
plot(checkvalve.t,checkvalve.Work_in_overtime,checkvalve.t,checkvalve.Work_out_overtime)
legend('Work in','Work out'), title(['Check Valve is ',num2str(checkvalve.e*100), ' % Efficient'])



%%
run('main_Scaled.m')
MYonepass
%% 
close all
figure(2)
labels = categorical({'EHA','HHEA','Check Valve'});
b = bar(labels,[e,Efficiency,checkvalve.e]);
ylabel('Net Efficiency')
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
%title(['Electric Generator Sizing ( Total Mean Power = ',num2str(-Total_Work_OverTime(end)/t(end)/1000),' kW)'])
ylim([0 4500])

figure(4)
labels = categorical({'EHA','HHEA','Check Valve'});
b = bar(labels,[EHA_Pump_size, 0; MP_size HECM_size; checkvalve.motorsize 0],'stacked');
ylabel('Size of hydraulic pump/motor (cc)')
legend('Main pump/motor','HECM pump/motor','Location','northwest')
xtips1 = b(2).XEndPoints;
ytips1 = b(2).YEndPoints;
labels1 = string(round(ytips1)) +  [' cc'];
labels1(2) = string(round(MP_size)) + [' cc +'] +  string(round(HECM_size)) + [' cc'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
%title('Hydraulic Pump/Motor Sizing')
ylim([0 1.25*(MP_size + HECM_size)])

figure(5)
plot(t, -F1, '-', t([1,end]),ones(2,1)*Frange1', '-', t, CPR_F1(d),'-o','linewidth',2); 
ylabel('Force (N)'); xlabel('Time (s)');grid
title('HHEA rail selections')

figure
P1 = plot(t,cumsum(-P_in)*dt/1000,t,cumsum(-P_out)*dt/1000,t,-Total_Work_OverTime/1000,checkvalve.t,checkvalve.Work_in_overtime/1000,checkvalve.t,checkvalve.Work_out_overtime/1000);
legend('Work in controlled','Electrical Energy out EHA','Electrical Energy out HHEA','Work in Check Valve','Electrical Energy Out Check Valve')
ylabel('Work (kJ)')
xlabel('Time (s)')
xlim([0 200])
%title('Energy Absorbtion')

figure
P1 = plot(t,cumsum(-P_in)*dt/1000,checkvalve.t,checkvalve.Work_in_overtime/1000);
legend('Work in PI control','Work in Check Valve')
ylabel('Work (kJ)')
xlabel('Time (s)')
xlim([0 200])
%title('Energy Absorbtion')

figure
P1 = plot(t,cumsum(-P_in)*dt,t,-Total_Work_OverTime,t,-MP_work,t,-HECM_work);
legend('Work in controlled','Total Electrical Energy out HHEA','Main Pump','HECM')
ylabel('Work (J)')
xlabel('Time (s)')
xlim([0 200])



PI_work_in = cumsum(-P_in)*dt;
EHA_work_out = cumsum(-P_out)*dt;
PI_avework_in = (PI_work_in(end) - PI_work_in(find(t==100)))/100; % W
EHA_avework_out = (EHA_work_out(end) - EHA_work_out(find(t==100)))/100; % W
HHEA_avework_out = (-Total_Work_OverTime(end) + Total_Work_OverTime(find(t==100)))/100; % W
checkvalve.avework_in = (checkvalve.Work_in_overtime(end) - checkvalve.Work_in_overtime(find(checkvalve.t==100)))/100; % W
checkvalve.avework_out = (checkvalve.Work_out_overtime(end) - checkvalve.Work_out_overtime(find(checkvalve.t==100)))/100; % W
MainPump_avework_out = (-MP_work(end) + MP_work(find(t==100)))/100; % W
HECM_avework_out = (-HECM_work(end) + HECM_work(find(t==100)))/100; % W


figure
labels = categorical({'Power in Check Valve','Check Valve out','EHA out','HHEA out','Power in PI'});
labels = reordercats(labels,{'Power in Check Valve','Check Valve out','EHA out','HHEA out','Power in PI'});
b = bar(labels,[checkvalve.avework_in/1000; checkvalve.avework_out/1000 ; EHA_avework_out/1000; HHEA_avework_out/1000; PI_avework_in/1000]);
ylabel('Average Power after ramp up (kW)')
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;
labels1 = string(round(ytips1)) +  [' kW'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
print('Figures/AvePower', '-dpng', '-r600')

figure
labels = categorical({'Check Valve','EHA','HHEA'});
labels = reordercats(labels,{'Check Valve','EHA','HHEA'});
b = bar(labels,[checkvalve.avework_in/1000 , checkvalve.avework_out/1000 ; PI_avework_in/1000, EHA_avework_out/1000; PI_avework_in/1000, HHEA_avework_out/1000;]);
ylabel('Average Power after ramp up (kW)')
legend('Average Power in','Average Power out')
print('Figures/AvePower', '-dpng', '-r600')

figure
plot(t,-P_in/1000,checkvalve.t,checkvalve.V.*checkvalve.F/1000)
legend('Power in PI control','Power in Check Valve')
ylabel('Power (kW)')
xlabel('Time (s)')
xlim([150 164])

%%
close all
figure
plot(t,-P_in/1000,t,-MP_power/1000,t,-battery/1000,'k')
xlabel('Time (s)')
ylabel('Power (kW)')
xlim([150 164])
grid on

%%
figure;
plot(t, -F1, '-',t([1,end]),ones(2,1)*Frange1', '-','linewidth',2)
hold on
scatter( t, CPR_F1(d),'o','MarkerEdgeColor',[0.6350 0.0780 0.1840]);
hold off
ylabel('Force (N)'); grid, xticks([]);%xlabel('Time (s)');
xlim([175 200])

figure;
plot(t, -F1, '-',t([1,end]),ones(2,1)*Frange1', '-','linewidth',2)
hold on
scatter( t, CPR_F1(d),'o','MarkerEdgeColor',[0.6350 0.0780 0.1840]);
hold off
ylabel('Force (N)'); xlabel('Time (s)');grid
xlim([157 158])
ylim([-.3 -.1]*1e7)
return

figure
plot(t,-P_in/1000)
ylabel('Power (kW)')
xlabel('Time (s)')
xlim([150 164])

figure
yyaxis left
plot(t,-P_in/1000,t,-MP_power/1000,t,-battery/1000,'k')
xlabel('Time (s)')
ylabel('Power (kW)')
%ylim([-2e5 12e5])
yyaxis right
plot(t,v,'Color','#D95319')
legend('Power In','Power Out Main Pump','Power Out HECM','Velocity','location','northwest')
xlim([175 200])
ylim([-2 12])
ylabel('Velocity (m/s)')
grid on
print('Figures/HHEAPower4', '-dpng', '-r600')

figure
plot(t,-P_in/1000,t,-MP_power/1000,t,-battery/1000,'k')
xlabel('Time (s)')
ylabel('Power (kW)')
legend('Power In','Power Out Main Pump','Power Out HECM','location','northwest')
xlim([150 164])
grid on


ActualMPLoss = sum(ActualMPLoss);
HECMLosses = sum(HECMLosses(d_ind))*dt;
TotalSwitchingLosses = TotalSwitchingLosses;
Totallosses = Totallosses;
figure
labels = categorical({'Main Pump Losses','HECM Losses','Switching Losses'});
b = bar(labels,[ActualMPLoss/1000, HECMLosses/1000, TotalSwitchingLosses/1000]);
ylabel('Energy Loss (kJ)')
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;
labels1 = string(round(ytips1)) +  [' kJ'];
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
print('Figures/LossContributions', '-dpng', '-r600')
%% EHA instantaneous efficiency

Hydraulic_in = (-P_in>=0).*(-P_in);
Hydraulic_out = -(-P_in<=0).*(-P_in);
Mechanical_in = -(-P_out<=0).*(-P_out);
Mechanical_out = (-P_out>=0).*(-P_out);
inst_e = (Hydraulic_out-Mechanical_out)./(Mechanical_in-Hydraulic_in);
%figure, plot(t,inst_e), ylabel('Instantaneous efficiency'), xlabel('Time (s)')


figure, plot(t,Hydraulic_in,t,-Hydraulic_out,t,Mechanical_in,t,-Mechanical_out), legend('Hydraulic in','Hydraulic out','Mechanical in','Mechanical out')

ave_Hydraulic_in = sum((-P_in>=0).*(-P_in))*dt;
ave_Hydraulic_out = sum(-(-P_in<=0).*(-P_in))*dt;
ave_Mechanical_in = sum(-(-P_out<=0).*(-P_out))*dt;
ave_Mechanical_out = sum((-P_out>=0).*(-P_out))*dt;
ave_inst_e = (ave_Hydraulic_out-ave_Mechanical_out)/(ave_Mechanical_in-ave_Hydraulic_in)
ave_inst_e = (ave_Hydraulic_out+ave_Mechanical_out)/(ave_Mechanical_in+ave_Hydraulic_in)


figure
plot(t,fracDisp)
ylabel('Fractional Displacement')
xlabel('Time (s)')
grid on

figure
plot(t,-P_in/1000,t,-P_out/1000)
xlabel('Time (s)')
ylabel('Power (kW)')
legend('Power In','Power Out','location','northwest')
xlim([150 164])
grid on
%% HHEA
ave_Hydraulic_in_HHEA = sum((-P_in>=0).*(-P_in))*dt;
ave_Hydraulic_out_HHEA = sum(-(-P_in<=0).*(-P_in))*dt;
ave_Mechanical_in_HHEA = sum(-(-battery<=0).*(-battery))*dt;
ave_Mechanical_out_HHEA = sum(-MP_power)*dt + sum((-battery>=0).*(-battery));
ave_inst_e_HHEA = (ave_Hydraulic_out_HHEA-ave_Mechanical_out_HHEA)/(ave_Mechanical_in_HHEA-ave_Hydraulic_in_HHEA)
ave_inst_e_HHEA = (ave_Hydraulic_out_HHEA+ave_Mechanical_out_HHEA)/(ave_Mechanical_in_HHEA+ave_Hydraulic_in_HHEA)

%% Rail Flow
close all
finalflow4 = sum(QR{4}(d_ind)*dt);
MPstaringt = 50; % seconds (must be a multiple of dt, otherwise the MP_work line of code needs to be changed
MP_flow = [zeros(find(t==MPstaringt),1); cumsum(finalflow4/(t(end)-MPstaringt)*ones(length(t)-find(t==MPstaringt),1))*dt];
    
figure
grid on
plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind)])*dt,t,MP_flow);
%legend([num2str(PR(1)*1e-6) ' MPa Pressure Rail'],[num2str(PR(2)*1e-6) ' MPa Pressure Rail'],[num2str(PR(3)*1e-6) ' MPa Pressure Rail'],[num2str(PR(4)*1e-6) ' MPa Pressure Rail'],'Main pump flow','location','northwest');
ylabel('Cumulative Rail Flow ($m^3$)','interpreter','latex');
xlabel('Time  ($s$)','interpreter','latex');
x = [.4 .5]; y = [.8 .7]; annotation('textarrow',x,y,'String',[num2str(PR(1)*1e-6,3) ' MPa Pressure Rail']);
x = [.55 .46]; y = [.61 .54]; annotation('textarrow',x,y,'String',[num2str(PR(2)*1e-6,3) ' MPa Pressure Rail']);
x = [.55 .5]; y = [.45 .495]; annotation('textarrow',x,y,'String',[num2str(PR(3)*1e-6,3) ' MPa Pressure Rail']);
x = [.55 .48]; y = [.2 .39]; annotation('textarrow',x,y,'String',[num2str(PR(4)*1e-6,3) ' MPa Pressure Rail']);
x = [.41 .45]; y = [.38 .435]; annotation('textarrow',x,y,'String','Flow from Main Pump');
grid on

    
figure
grid on
plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind)])*dt,t,MP_flow);
%legend([num2str(PR(1)*1e-6) ' MPa Pressure Rail'],[num2str(PR(2)*1e-6) ' MPa Pressure Rail'],[num2str(PR(3)*1e-6) ' MPa Pressure Rail'],[num2str(PR(4)*1e-6) ' MPa Pressure Rail'],'Main pump flow','location','northwest');
ylabel('$m^3$','interpreter','latex');
xlabel('Time ($s$)','interpreter','latex');
xlim([110 125])
ylim([-2.5 -1.5])
grid on


figure
subplot(3,3,[4 5 7 8])
grid on
plot(t, cumsum([QR{1}(d_ind),QR{2}(d_ind),QR{3}(d_ind),QR{4}(d_ind)])*dt,t,MP_flow);
grid on
%yticks(-5:2.5:5)
%legend([num2str(PR(1)*1e-6) ' MPa Pressure Rail'],[num2str(PR(2)*1e-6) ' MPa Pressure Rail'],[num2str(PR(3)*1e-6) ' MPa Pressure Rail'],[num2str(PR(4)*1e-6) ' MPa Pressure Rail'],'Main pump flow','location','northwest');
ax = gca;
area = [80 -1.8 100 -.8];
inlarge = subplot(3,3,3);
panpos = inlarge.Position;
delete(inlarge);
inlarge = zoomin(ax,area,panpos);
grid on
title(inlarge,'Zoom in')

function pan = zoomin(ax,areaToMagnify,panPosition)
% AX is a handle to the axes to magnify
% AREATOMAGNIFY is the area to magnify, given by a 4-element vector that defines the
%      lower-left and upper-right corners of a rectangle [x1 y1 x2 y2]
% PANPOSTION is the position of the magnifying pan in the figure, defined by
%        the normalized units of the figure [x y w h]
%

fig = ax.Parent;
pan = copyobj(ax,fig);
pan.Position = panPosition;
pan.XLim = areaToMagnify([1 3]);
pan.YLim = areaToMagnify([2 4]);
pan.XTick = [];
pan.YTick = [];
rectangle(ax,'Position',...
    [areaToMagnify(1:2) areaToMagnify(3:4)-areaToMagnify(1:2)])
xy = ax2annot(ax,areaToMagnify([1 4;3 2]));
annotation(fig,'line',[xy(1,1) panPosition(1)],...
    [xy(1,2) panPosition(2)+panPosition(4)],'Color','k')
annotation(fig,'line',[xy(2,1) panPosition(1)+panPosition(3)],...
    [xy(2,2) panPosition(2)],'Color','k')
end

function anxy = ax2annot(ax,xy)
% This function converts the axis unites to the figure normalized unites
% AX is a handle to the figure
% XY is a n-by-2 matrix, where the first column is the x values and the
% second is the y values
% ANXY is a matrix in the same size of XY, but with all the values
% converted to normalized units

pos = ax.Position;
%   white area * ((value - axis min) / axis length)   + gray area
normx = pos(3)*((xy(:,1)-ax.XLim(1))./range(ax.XLim))+ pos(1);
normy = pos(4)*((xy(:,2)-ax.YLim(1))./range(ax.YLim))+ pos(2);
anxy = [normx normy];
end
