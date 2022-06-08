clear, close all
kp_0 = 2.9897e7;
ki_0 = 3.2542e7;
x0 = [kp_0;ki_0];

tStart = tic;
return
options = optimset('PlotFcns',@optimplotfval);
%[x,fval] = fminsearch(@fun,x0,options)

% nvars = 2; options = optimoptions('ga','PlotFcn', @gaplotbestf);
% [x,fval] = ga(@fun,nvars,[],[],[],[],zeros(nvars,1),[],[],options);

TimeElapsed = toc(tStart)  

%%
N = 25;
k_vals = linspace(5e6,5e7,N);

KP = NaN(N,N); KI = KP; J = KP;
for i = 1:N
    (i-1)*100/N
    for j = 1:N
        KP(i,j) = k_vals(i);
        KI(i,j) = k_vals(j);
        J(i,j) = fun([k_vals(i) ; k_vals(j)]);
    end
end
figure, contour(KP,KI,J), xlabel('K_p'), ylabel('K_i')

%%
kp = 1.9e07;                             % PTO Damping Coeff [Nsm/rad]
ki = 3.3e7; 
fun([kp ki])

function J = fun(x)
kp = x(1);
ki = x(2);
%A = x(3);
[~]=evalc('wecSim');
EHA_losses;
%J = Work_Out;
J = Work_in;
end