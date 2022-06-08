Pmax_set = [28.2 28.5 28.8 29.1 29.4 29.7 30 30.6 30.9 31.5 31.8 32];

T1_set = [11 11 11 11 12 12 12 12 12 12 12 13];
T2_set = [20 20 20 20 20 21 21 21 21 22 22 22];
T3_set = [14 14 15 15 15 14 15 15 15 15 16 15];
T4_set = [19 20 20 20 20 20 21 21 21 22 22 22];

eff_nd = [104.67 104.81 104.95 104.91 104.98 105.04 105.43 105.31...
    105.22 105.67 105.81 105.83];
eff_gr = [70.22 70.21 70.37 70.36 70.47 70.72 70.93 70.92 70.94...
    71.53 71.8 71.73];
eff_tr = [87.10 87.09 87.17 87.13 87.2 87.16 87.33 87.29 87.24...
    87.37 87.45 87.4];

% figure(6)
% subplot(2,2,1);
% plot(Pmax_set,T1_set)
% xlabel('Maximum Pressure')
% ylabel('Maximum Boom Torque')
% subplot(2,2,2);
% plot(Pmax_set,T2_set)
% xlabel('Maximum Pressure')
% ylabel('Maximum Arm Torque')
% subplot(2,2,3);
% plot(Pmax_set,T3_set)
% xlabel('Maximum Pressure')
% ylabel('Maximum Bucket Torque')
% subplot(2,2,4);
% plot(Pmax_set,T4_set)
% xlabel('Maximum Pressure')
% ylabel('Maximum Swing Torque')

figure(7)
plot(Pmax_set,(T1_set+T2_set+T3_set+T4_set))
xlabel('Maximum Pressure')
ylabel('Addition of Torque Limits')

figure(8)
subplot(2,2,1);
plot(Pmax_set,eff_nd)
xlabel('Maximum Pressure')
ylabel('Efficiencies - Ninety Degree')

subplot(2,2,2);
plot(Pmax_set,eff_gr)
xlabel('Maximum Pressure')
ylabel('Efficiencies - Grading')

subplot(2,2,3);
plot(Pmax_set,eff_tr)
xlabel('Maximum Pressure')
ylabel('Efficiencies - Trenching')