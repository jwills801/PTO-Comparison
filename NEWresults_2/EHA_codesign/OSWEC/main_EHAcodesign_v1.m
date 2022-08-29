Np = 10; kp_vals = linspace(.5e7,5e7,Np);
Ni = 10; ki_vals = linspace(.5e7,5e7,Ni);

KP = NaN(Np,Ni); KI = KP; J = KP;
outertime = tic;
for i = 1:Np
    for j = 1:Ni
        KP(i,j) = kp_vals(i);
        KI(i,j) = ki_vals(j);
        J(i,j) = fun([kp_vals(i) ; ki_vals(j)]);
        disp([num2str((i-1)*Np+j) ' of ' num2str(Np*Ni) ' simulations complete'])
    end
end
toc(outertime)
figure, contour(KP,KI,J,25), xlabel('K_p'), ylabel('K_i')

KP_ = KP(:); KI_ = KI(:); J_ = J(:);
[a,b] = min(J_);
disp('Best case was:')
disp(['              ' num2str(-a/1e6) ' MJ'])
disp(['    with Kp = ' num2str(KP_(b)/1e7) 'e7'])
disp(['     and Ki = ' num2str(KI_(b)/1e7) 'e7'])

%figure, plot(KP_,KI_,'*')


return
%%
tic
kp = 1.9e07;                             % PTO Damping Coeff [Nsm/rad]
ki = 3.3e7; 
fun([kp ki])
toc

function J = fun(x)
kp = x(1);
ki = x(2);
%A = x(3);
[~]=evalc('wecSim');
EHA_losses;
%J = Work_Out;
J = Work_In;
end