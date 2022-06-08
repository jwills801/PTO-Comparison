
%Example of user input MATLAB file for post processing
close all

%Plot waves
waves.plotEta(simu.rampTime);
try 
    waves.plotSpectrum();
catch
end

% Plot RY forces for body 1
plotForces(output,1,5)

%Plot RY response for body 1
output.plotResponse(1,5);

% Plot x forces for body 2
plotForces(output,2,1)

figure(), plot(output.ptos.time,output.ptos.powerInternalMechanics), ylabel('Power (W)'), title('powerInternalMechanics')
figure(), plot(output.ptos.time,cumsum(output.ptos.powerInternalMechanics)*simu.dt), ylabel('Energy (J)'), title('Energy InternalMechanics')

ramp_ind = find(simu.time == simu.rampTime);
Ave_Absorbed_Power = -mean(output.ptos.powerInternalMechanics(ramp_ind:end,5));
CWR = Ave_Absorbed_Power/waves.Pw/18

