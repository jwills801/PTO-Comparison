ramp_ind = find(simu.time == 50);
Ave_Absorbed_Power = -mean(output.ptos.powerInternalMechanics(ramp_ind:end,5));
CWR(imcr) = Ave_Absorbed_Power/waves.Pw/18;

periods(imcr) = waves.T;


%% Code to run after the mcr run to get CWR plot overlaid with the irregular wave spectrum
load PM_T8_H25
figure
yyaxis left
plot(periods,CWR), ylabel('Capture Width Ratio'), xlabel('Period (s)')
yyaxis right
plot(2*pi./w,S) % w and S are variables saved in matlab file containing irregular spectrum
xlim([1 13])
grid on
ylabel('Spectrum (m^2 s/rad)')