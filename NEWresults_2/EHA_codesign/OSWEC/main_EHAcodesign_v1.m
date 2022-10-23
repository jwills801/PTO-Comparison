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


%% Find best PI gains for multiple frequencies, then try to make a transfer funtion to approximate them all at their frequencies
clear, close all
Np = 7; kp_vals = linspace(.25e7,10e7,Np);
Ni = 7; ki_vals = linspace(.5e7,7e7,Ni);
Nw = 5; w_vals = linspace(.55,1.5,Nw);
% reference values from PM spectrum with H = 2.5 and T = 8
w = [.5260 .5522 .6554 .7607  .8671  .9720  1.0774 1.1831 1.2903 1.3958 1.4994 3.9265];
A = [.0739 .1737 .9332 1.4098 1.3072 1.0056 0.7192 0.5030 .35    .2475  0.1785 0.0016];

OUTERtime = tic;
for k = 1:Nw
    T = 2*pi/w_vals(k);
    H = 2*interp1(w,A,w_vals(k));
    KP = NaN(Np,Ni); KI = KP; J = KP;
    outertime = tic;
    for i = 1:Np
        for j = 1:Ni
            KP(i,j) = kp_vals(i);
            KI(i,j) = ki_vals(j);
            J(i,j) = fun([kp_vals(i) ; ki_vals(j) ; T ; H]);
            disp([num2str((i-1)*Np+j) ' of ' num2str(Np*Ni) ' simulations complete (Frequency ',num2str(k),' of ',num2str(Nw),')'])
        end
    end
    toc(outertime)
    figure, contour(KP,KI,J,25), xlabel('K_p'), ylabel('K_i')

    KP_ = KP(:); KI_ = KI(:); J_ = J(:);
    [a,b] = min(J_);
    disp(['Best case for regualar waves with period T = ',num2str(round(T,2)),' seconds:'])
    disp(['              ' num2str(-a/1e6) ' MJ'])
    disp(['    with Kp = ' num2str(KP_(b)/1e7) 'e7'])
    disp(['     and Ki = ' num2str(KI_(b)/1e7) 'e7'])
    Best_Kps(k) = KP_(b);
    Best_Kis(k) = KI_(b);
    Tf = KP_(b) + KI_(b)/sqrt(-1)/w_vals(k);
    Phase(k) = angle(Tf)*180/pi; % degrees
    Mag_dB(k) = 20*log10(abs(Tf)); % dB
end
figure
subplot(211), semilogx(w_vals,Mag_dB,'*'), ylabel('Magnitude [dB]'), grid, xlim([.1 10])
subplot(212), semilogx(w_vals,Phase,'*'), xlabel('Frequency [rad/s]'), ylabel('Phase [Degrees]'), grid, xlim([.1 10])

toc(OUTERtime)

% RESULTS
% Took 3 hours
%w_vals = [0.5500    0.6556    0.7611    0.8667    0.9722    1.0778    1.1833    1.2889    1.3944    1.5000];
%Phase = [     -69.9795  -71.8520  -66.0981  -58.0237  -52.4804  -41.7866  -31.2409  -21.9640  -17.2459  -12.2307]; % Degrees
%Mag_dB = [    148.8322  149.6534  152.1592  155.1668  157.2325  157.8501  157.6428  157.8183  156.6817  155.5000]; % dB
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
T = x(3);
H = x(4);
[~]=evalc('wecSim');
EHA_losses;
%J = Work_Out;
J = Work_In;
end