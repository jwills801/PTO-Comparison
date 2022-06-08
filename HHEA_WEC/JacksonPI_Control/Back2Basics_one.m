% Implement a PI controller on the simple system
% Parameters are taken from WECSIM

% M vdot + B v + k int(v dt) = Fe - Fpto
load('Fexc.mat')
M = 8.5677e6;
B = 2.0776e7;
K = 1.9533e5;
dt = t(2)-t(1);
w = 2*pi/8; % rad/s

% simulate without control
Fpto = zeros(size(t));
v = 0;
x = 0;
for i = 1:length(t)-1
    vdot = 1/M*(Fexc_Pitch(i)- Fpto(i) - B*v(i) - K*x(i));
    v(i+1) = v(i) + vdot*dt;
    x(i+1) = x(i) + v(i)*dt;
end
figure(), plot(t,Fpto,t,v), legend('PTO force','vel'), xlabel('Time (s)'), ylabel('N or m/s'), title(['No PTO force => No Energy absorbed'])



% Simulate w damping control
v = 0;
x = 0;
for i = 1:length(t)-1
    Fpto(i) = -B*v(i);
    vdot = 1/M*(Fexc_Pitch(i)- Fpto(i) - B*v(i) - K*x(i));
    v(i+1) = v(i) + vdot*dt;
    x(i+1) = x(i) + v(i)*dt;
end
figure(), plot(t,Fpto/1000/1000,t,v*10)
legend('PTO force','vel'), xlabel('Time (s)'), ylabel('MN or cm/s'),title(['Damping Control => Energy absorbed: ',num2str(sum(-Fpto.*v')/1e9),' GJ'])
%xlim([50 100])



% Simulate w PI control
s = tf('s');
C = B + (w^2*M - K)/s;
C_discrete = c2d(C,dt)
[z, g] = zero(C_discrete);
p = pole(C_discrete);
v = 0;
x = 0;
for i = 1:length(t)-1   
    vdot = 1/M*(Fexc_Pitch(i)- Fpto(i) - B*v(i) - K*x(i));
    v(i+1) = v(i) + vdot*dt;
    x(i+1) = x(i) + v(i)*dt;
    Fpto(i+1) = g*(v(i) - z*v(i+1)) + p*Fpto(i);
end
figure(), plot(t,Fpto/1000/1000,t,v*10)
legend('PTO force','vel'), xlabel('Time (s)'), ylabel('MN or cm/s'),title(['PI Control => Energy absorbed: ',num2str(sum(-Fpto.*v')/1e9),' GJ'])
%xlim([50 100])

figure(), plot(t,cumsum(-Fpto.*v')), ylabel('Energy Absorbed'), xlabel('Time (s)'), title('PI control')






%% Hello
options = optimoptions(@ga,'PopulationSize',30,'MaxGenerations',5);
[x,fval] = ga(@(K)pitest(K),2,[],[],[],[],[],[],[],options);
x
1/fval



function J = pitest(K)
Kp = K(1) * 1e6;
Ki = K(2) * 1e6;
set_param('Back2Basics_two/PID Controller','P',num2str(Kp))
set_param('Back2Basics_two/PID Controller','I',num2str(Ki))
out = sim('Back2Basics_two');
dt = diff(out.tout);
[Kp Ki]
energy = sum(out.power)
plot(out.tout,cumsum(out.power)./[dt; dt(end)])
drawnow
J = 1/energy;
if energy < 0
    J = inf;
end
end

