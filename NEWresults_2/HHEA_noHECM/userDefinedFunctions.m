
%Plot waves
waves.plotEta(simu.rampTime);
try 
    waves.plotSpectrum();
catch
end



v = squeeze(myoutput.signals.values(11,:,:));
F = squeeze(myoutput.signals.values(17,:,:));
time = myoutput.time;
figure(), plot(time,v,time,F/1e6), ylabel('m/s or MN'), title('Linear PTO Force and velocity')

%%
figure(), plot(time,-cumsum(v.*F)*simu.dt/1e6), ylabel('Mechanical Energy (MJ)'), title('Absorbed energy at actuator'), xlabel('Time (s)')


figure, yyaxis left, plot(time,F/1e6), ylabel('Force [MN]'), xlabel('Time [s]'), ylim([-ceil(max(abs(F/1e6))) ceil(max(abs(F/1e6)))]), grid on
yyaxis right, plot(time,v), ylabel('Velocity [m/s]'), xlim([100 125]), ylim([-ceil(max(abs(v))) ceil(max(abs(v)))]), grid on
