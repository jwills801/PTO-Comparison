return
%Example of user input MATLAB file for post processing
%return
close all

%Plot waves
waves.plotEta(simu.rampTime);
try 
    waves.plotSpectrum();
catch
end




return

% Plot RY forces for body 1
plotForces(output,1,5)

%Plot RY response for body 1
output.plotResponse(1,5);

% Plot x forces for body 2
plotForces(output,2,1)

figure(), plot(output.ptos.time,output.ptos.powerInternalMechanics), ylabel('Power (W)'), title('powerInternalMechanics')
figure(), plot(output.ptos.time,cumsum(output.ptos.powerInternalMechanics)*simu.dt), ylabel('Energy (W)'), title('Energy InternalMechanics')


translational_vel = myoutput.signals.values(:,11);
translational_pto_force = myoutput.signals.values(:,17);
time = myoutput.time;
figure(), plot(time,translational_vel,time,translational_pto_force/1e6), ylabel('m/s or MN'), title('Linear PTO Force and velocity')
figure(), plot(time,cumsum(translational_vel.*translational_pto_force)*simu.dt), ylabel('J'), title('Absorbed energy at actuator')
