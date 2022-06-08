% %Example of user input MATLAB file for post processing
% close all
% 
% %Plot waves
% waves.plotEta(simu.rampTime);
% try 
%     waves.plotSpectrum();
% catch
% end
% 
% % Plot RY forces for body 1
% plotForces(output,1,5)
i = 5;
plot(output.bodies(1).time,-cumsum(output.bodies(1).velocity(:,i).*output.ptos.forceTotal(:,i)))
% 
% %Plot RY response for body 1
% output.plotResponse(1,5);
% 
% % Plot x forces for body 2
% plotForces(output,2,1)







% 5500