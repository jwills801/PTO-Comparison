% Calculate PI controller from data
wecSim;
l = 5.5; period = 8; freq = 1/period;
Fex_ = output.bodies(1).forceExcitation;
vel_ = output.bodies(1).velocity;
pos_ = output.bodies(1).position;
x = pos_(:,1); z = pos_(:,3)+10;
t = output.bodies(1).time;
dt = mean(diff(t));

Fex = Fex_(:,1).*pos_(:,3)-Fex_(:,3).*pos_(:,1) + Fex_(:,5);
w_check = (vel_(:,1).*z-vel_(:,3).*x)/l^2;
w = vel_(:,5);
plot(t,[w,w_check]);

%% Find Kp and KI from openloop response

[C,lags]=xcorr(w,-Fex,40);
[mC,ind]=max(C);
plot(lags*dt,C); grid

gain = (max(Fex)-min(Fex))/(max(w)-min(w));
delay = lags(ind)*dt;
plot(t+delay,-Fex, t,gain*w) 
phas = (delay/period*2*pi)

Kp = gain/sqrt(1+tan(phas)^2)
KI = Kp*(freq*2*pi)*tan(phas)
