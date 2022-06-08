wecSim;

X = output.bodies(1).position(1:end-1,5); 
T = simu.dt;                      % Sampling period       
Fs = 1/T;                         % Sampling frequency                    
L = length(simu.time(1:end-1));   % Length of signal

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure
yyaxis left
plot(2*pi*f,P1), xlim([0 3])
title('Single-Sided Amplitude Spectrum of Pitch Angle')
xlabel('Frequency (rad/s)')
ylabel('|P1(f)|')
yyaxis right
plot(waves.w,waves.S)
ylabel('Spectrum (m^2 s/rad)')