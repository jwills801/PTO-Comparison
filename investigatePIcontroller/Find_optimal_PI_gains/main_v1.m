clear, close all
kp_0 = 2.9897e7;
ki_0 = 3.2542e7;
x0 = [kp_0;ki_0];

tStart = tic;

options = optimset('PlotFcns',@optimplotfval);
[x,fval] = fminsearch(@fun,x0,options)


TimeElapsed = toc(tStart)  




function Work_in = fun(x)
kp = x(1);
ki = x(2);
%A = x(3);
[~]=evalc('wecSim');
ramp_ind = find(simu.time == simu.rampTime);
Work_in = sum(output.ptos.powerInternalMechanics(ramp_ind:end,5))*simu.dt;
end