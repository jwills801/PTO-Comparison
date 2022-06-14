
%Plot waves
waves.plotEta(simu.rampTime);
try 
    waves.plotSpectrum();
catch
end



translational_vel = squeeze(myoutput.signals.values(11,:,:));
translational_pto_force = squeeze(myoutput.signals.values(17,:,:));
time = myoutput.time;
figure(), plot(time,translational_vel,time,translational_pto_force/1e6), ylabel('m/s or MN'), title('Linear PTO Force and velocity')

%%
figure(), plot(time,-cumsum(translational_vel.*translational_pto_force)*simu.dt/1e6), ylabel('Mechanical Energy (MJ)'), title('Absorbed energy at actuator'), xlabel('Time (s)')
