clear, close all
load noPTO_nov3_2021.mat

if simu.morisonElement == 0 & simu.yawNonLin == 0 & simu.b2b == 0 & simu.nlHydro == 0
else
    disp('The equation of motion are different from what are assumed here, results will not be good.')
end

a =5; b = 3.9;

%% Radiation Damping
F_radDamp_calc = (body(1).hydroForce.fDamping * output.bodies(1).velocity')';

% figure, plot(output.bodies(1).time,output.bodies(1).forceRadiationDamping)
% figure, plot(output.bodies(1).time,F_radDamp_calc)
% figure, plot(output.bodies(1).time,F_radDamp_calc-output.bodies(1).forceRadiationDamping)


syms t
syms x(t)

XDOT = [a*diff(x,t)*cos(x);0;a*diff(x,t)*sin(x);0;diff(x,t);0];
rad_func = [a*cos(x) 0 -a*sin(x) 0 1 0]*body(1).hydroForce.fDamping*XDOT



diff(x(t), t)*((9519113892690845*sin(x(t)))/36028797018963968 - (37618636679434265*cos(x(t)))/17179869184 + 3355659414200545/17179869184) + ... 
    5*cos(x(t))*diff(x(t), t)*((5311903631544805*cos(x(t)))/1073741824 + (20152136796327855*sin(x(t)))/72057594037927936 - 7581376454787105/17179869184) ... 
    - 5*sin(x(t))*diff(x(t), t)*((8675524737733615*sin(x(t)))/137438953472 - (30514789644207155*cos(x(t)))/2251799813685248 + 1704296393812919/1125899906842624)


-37618636679434265/17179869184 + 3355659414200545/17179869184 + 5*5311903631544805/1073741824 - 5*7581376454787105/17179869184





%diff(x(t), t)*(((9519113892690845*sin(x(t)))/36028797018963968 - (37618636679434265*cos(x(t)))/17179869184 + 3355659414200545/17179869184) + ...
%    5*cos(x(t))*((5311903631544805*cos(x(t)))/1073741824 + (20152136796327855*sin(x(t)))/72057594037927936 - 7581376454787105/17179869184) - ...
%    5*sin(x(t))*((8675524737733615*sin(x(t)))/137438953472 - (30514789644207155*cos(x(t)))/2251799813685248 + 1704296393812919/1125899906842624));

eval( subs( (((9519113892690845*sin(x(t)))/36028797018963968 - (37618636679434265*cos(x(t)))/17179869184 + 3355659414200545/17179869184) + ...
    5*cos(x(t))*((5311903631544805*cos(x(t)))/1073741824 + (20152136796327855*sin(x(t)))/72057594037927936 - 7581376454787105/17179869184) - ...
    5*sin(x(t))*((8675524737733615*sin(x(t)))/137438953472 - (30514789644207155*cos(x(t)))/2251799813685248 + 1704296393812919/1125899906842624)),x,0))

% THERE IS NO LINEAR DAMPING OR MORISON ELEMENT


%% Quadratic damping
XDOT2 = [a*diff(x,t)*cos(x);0;a*diff(x,t)*sin(x);0;diff(x,t);0].^2;
quad_damp_func = [a*cos(x) 0 -a*sin(x) 0 1 0]*body(1).hydroForce.visDrag*XDOT2


(7034159837642489*diff(x(t), t)^2)/268435456 + (755970078670847875*cos(x(t))^3*diff(x(t), t)^2)/68719476736 - (3736125*sin(x(t))^3*diff(x(t), t)^2)/2

% IF this function is linearlized about x=0 and diff(x,t)=0, the result is 0


% THE Kp I get from this is 2.0535e7 but the real result is:

    % 2.0188e7





%% GET Ki term

% Stiffness
X = [a*sin(x);0;-a-b+a*cos(x);0;x;0];
stiff_func = [a*cos(x) 0 -a*sin(x) 0 1 0]*body(1).hydroForce.linearHydroRestCoef*X


x(t)*((5808321602479115*sin(x(t)))/18014398509481984 - 8512151767901799/4294967296)+ ...
    (89/10 - 5*cos(x(t)))*((27985156497211395*sin(x(t)))/17179869184 + 1161664320495823/18014398509481984)
 
 
%linearized:
-8512151767901799/4294967296 - 5*27985156497211395/17179869184 + 27985156497211395/17179869184*89/10

4.3710e+06

% Plus buoyancy:
(1271000-1744000)*(-a*sin(x))

Stiffness = 2365000 + 4.3710e+06


%% Added Mass and mass

XDDOT = [a*diff(x,t,t)*cos(x) - a*diff(x,t)^2*sin(x);0;-a*diff(x,t,t)*sin(x) - a*diff(x,t)^2*cos(x);0;diff(x,t,t);0];
M_matrix = diag([body(1).mass,body(1).mass,body(1).mass,[body(1).momOfInertia]]);
AddedMass_func = [a*cos(x) 0 -a*sin(x) 0 1 0]*body(1).hydroForce.fAddedMass*XDDOT + [a*cos(x) 0 -a*sin(x) 0 1 0]*M_matrix*XDDOT






(5*cos(x(t))*diff(x(t), t)^2 + 5*sin(x(t))*diff(x(t), t, t))*((14170884365126775*sin(x(t)))/34359738368 - (20052052723586205*cos(x(t)))/1125899906842624 + 3672588802909549/562949953421312) + ... 
    (5*sin(x(t))*diff(x(t), t)^2 - 5*cos(x(t))*diff(x(t), t, t))*((3357289052768695*sin(x(t)))/9007199254740992 - (32621689341109075*cos(x(t)))/2147483648 + 2931553420784817/2147483648) + ... 
    diff(x(t), t, t)*((24827802970182165*sin(x(t)))/18014398509481984 - (29170052075018435*cos(x(t)))/4294967296 + 7213064817628235/1073741824) - ... 
    635000*cos(x(t))*(5*sin(x(t))*diff(x(t), t)^2 - 5*cos(x(t))*diff(x(t), t, t)) + 635000*sin(x(t))*(5*cos(x(t))*diff(x(t), t)^2 + 5*sin(x(t))*diff(x(t), t, t)) + 1850000*diff(x(t), t, t)






5*32621689341109075/2147483648 - 5*2931553420784817/2147483648 -29170052075018435/4294967296 + 7213064817628235/1073741824 + 635000*5 + 1850000

7.4079e+07

%% Ki term
0.7854^2*7.4079e+07 - 6736000

3.8960e+07

% We want it to be 3.8278e7




return

%% Added mass
F_addedMass_calc = (body(1).hydroForce.fAddedMass * output.bodies(1).acceleration')';

% figure, plot(output.bodies(1).time,output.bodies(1).forceAddedMass)
% figure, plot(output.bodies(1).time,F_addedMass_calc)
% figure, plot(output.bodies(1).time,output.bodies(1).forceAddedMass-F_addedMass_calc)

%% Does decomposing the total force, result in an acceleration that matches the real results?
realF_total = output.bodies(1).forceTotal + output.ptos.forceConstraint;
a =5; b = 3.9;
% F_total_1DOF = output.bodies(1).forceTotal(:,5) + output.bodies(1).forceTotal(:,1)*a.*cos(output.bodies(1).velocity(:,5)) - output.bodies(1).forceTotal(:,3)*a.*sin(output.bodies(1).velocity(:,5));
F_total_1DOF = realF_total(:,5) + realF_total(:,1)*a.*cos(output.bodies(1).position(:,5)) - realF_total(:,3)*a.*sin(output.bodies(1).position(:,5));

F_total_1DOF_FuncOfTheta = body(1).mass*a^2*cos(output.bodies(1).position(:,5)).*(output.bodies(1).acceleration(:,5).*cos(output.bodies(1).position(:,5)) - (output.bodies(1).velocity(:,5)).^2.*sin(output.bodies(1).position(:,5))) + ...
    body(1).mass*a^2*sin(output.bodies(1).position(:,5)).*(output.bodies(1).acceleration(:,5).*sin(output.bodies(1).position(:,5)) + (output.bodies(1).velocity(:,5)).^2.*cos(output.bodies(1).position(:,5))) + ...
    body(1).momOfInertia(2)*output.bodies(1).acceleration(:,5);

figure, plot(output.bodies(1).time,F_total_1DOF)

figure, plot(output.bodies(1).time,F_total_1DOF_FuncOfTheta)

M_matrix = diag([body(1).mass,body(1).mass,body(1).mass,[body(1).momOfInertia]]);
figure, plot(output.bodies(1).time,M_matrix*output.bodies(1).acceleration'), title('Mass Times Acceleration'), ylabel('Force (N)'), xlabel('Time (s)')

figure, plot(output.bodies(1).time,realF_total), title('Sum of Forces'), ylabel('Force (N)'), xlabel('Time (s)')
figure, plot(output.bodies(1).time,F_total_1DOF_FuncOfTheta,output.bodies(1).time,F_total_1DOF), legend('Ftotal Function of theta','Ftotal 1DOF')

%%
x = a*sin(output.bodies(1).position(:,5));
z =-a-b+ a*cos(output.bodies(1).position(:,5));
t= output.ptos.time;
figure, plot(t,x,t,output.bodies(1).position(:,1))
figure, plot(t,z,t,output.bodies(1).position(:,3))
figure, plot(t,output.bodies(1).position(:,5))


%%
w = 2*pi/waves.T; % rad/s
M_total_matrix = M_matrix + body(1).hydroForce.fAddedMass;
M_total_scalar = M_total_matrix(1,1)*a^2 + 2*M_total_matrix(1,5)*a+M_total_matrix(5,5);

S_matrix = body(1).hydroForce.linearHydroRestCoef;
S_scalar = S_matrix(1,1)*a^2 + S_matrix(1,5)*a - S_matrix(3,3)*a^2 + S_matrix(3,3)*a*(a+b) + S_matrix(5,1)*a + S_matrix(5,5);

Ki = -w^2*M_total_scalar + S_scalar




%pto(1).k = 3.8278e7;                                 % PTO Stiffness Coeff [Nm/rad]
%pto(1).c = 2.0188e7;                             % PTO Damping Coeff [Nsm/rad]





