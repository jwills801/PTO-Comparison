clear
load hydrodata.mat
n = 6;
M = Mass(1:n,1:n) + addedmass(1:n,1:n);
C = LinearDamping(1:n,1:n) + radiationDamping(1:n,1:n);
K = hydrostaticstiffness(1:n,1:n);
Fexc = Fexc(:,1:n);
rho = 1025;
g = 9.81;


syms theta2 theta2dot theta2ddot
a =5; b = 3.9;
x_m = [a*sin(theta2);0;-a-b+a*cos(theta2);0;theta2;0];
x_mdot = [a*theta2dot*cos(theta2);0;-a*theta2dot*sin(theta2);0;theta2dot;0];
x_mddot = [-a*theta2dot^2*sin(theta2)+a*theta2ddot*cos(theta2);0;-a*theta2dot^2*cos(theta2)-a*theta2ddot*sin(theta2);0;theta2ddot;0];



Fexc5 = Fexc(:,5) + Fexc(:,1)*a.*cos(x(:,5)) - Fexc(:,3)*a.*sin(x(:,5));


Fk = K(5,5).*x(:,5) - K(3,3).*sin(x(:,5)).*(-a-b+a.*cos(x(:,5)));


Fc = C(1,4).*a.*vel(:,5).*cos(x(:,5)) + C(5,5).*vel(:,5) + C(1,1).*vel(:,5).*a^2.*(cos(x(:,5))).^2 + C(1,5).*vel(:,5).*a.*cos(x(:,5)) + C(3,3).*vel(:,5).*a^2.*(sin(x(:,5))).^2;
Fm = M(5,1)*( a.*acc(:,5).*cos(x(:,5)) - a.*vel(:,5).^2.*sin(x(:,5)) ) + M(5,5).*acc(:,5) + M(1,1).*a.*cos(x(:,5)).*(   a.*acc(:,5).*cos(x(:,5)) - a.*vel(:,5).^2.*sin(x(:,5))  ) + M(1,5).*acc(:,5).*a.*cos(x(:,5)) + M(3,3).*a.*sin(x(:,5)) .* ( a.*vel(:,5).^2.*cos(x(:,5)) + a.*acc(:,5).*sin(x(:,5))  );

LHS = Fk + Fc + Fm;

figure, plot(t,LHS,t,Fexc5), legend('f($\theta$)','F\_{exc}','Interpreter','Latex'), ylabel('Force (N)'), xlabel('Time(s)')
%figure, plot(t,LHS-Fexc5), ylabel('Force (N)'), xlabel('Time(s)')

%% Lets try each force individually:

%% ADDED MASS
F_addedmass = output.bodies(1).forceAddedMass;

% In our 1 DOF
F_addedmass5 = F_addedmass(:,5) + F_addedmass(:,1)*a.*cos(x(:,5)) - F_addedmass(:,3)*a.*sin(x(:,5));

% compare to M_add * X
Madd = addedmass(1:n,1:n);
F_MaddX = Madd(5,1)*( a.*acc(:,5).*cos(x(:,5)) - a.*vel(:,5).^2.*sin(x(:,5)) ) + Madd(5,5).*acc(:,5) + Madd(1,1).*a.*cos(x(:,5)).*(   a.*acc(:,5).*cos(x(:,5)) - a.*vel(:,5).^2.*sin(x(:,5))  ) + Madd(1,5).*acc(:,5).*a.*cos(x(:,5)) + Madd(3,3).*a.*sin(x(:,5)) .* ( a.*vel(:,5).^2.*cos(x(:,5)) + a.*acc(:,5).*sin(x(:,5))  );

figure, plot(t,F_addedmass5,t,F_MaddX), legend('From Output','My model')
figure, plot(t,F_addedmass5)
figure, plot(t,F_MaddX)

% Checks out!

%% Linear Damping
F_linearDamping = output.bodies(1).forceLinearDamping;

% in out 1 DOF
F_linearDamping5 = F_linearDamping(:,5) + F_linearDamping(:,1)*a.*cos(x(:,5)) - F_linearDamping(:,3)*a.*sin(x(:,5));

% Compare to C * Xdot
C = LinearDamping(1:n,1:n);
F_C_lin = C(1,4).*a.*vel(:,5).*cos(x(:,5)) + C(5,5).*vel(:,5) + C(1,1).*vel(:,5).*a^2.*(cos(x(:,5))).^2 + C(1,5).*vel(:,5).*a.*cos(x(:,5)) + C(3,3).*vel(:,5).*a^2.*(sin(x(:,5))).^2;


figure, plot(t,F_linearDamping5,t,F_C_lin), legend('From Output','My model')
figure, plot(t,F_linearDamping5)
figure, plot(t,F_C_lin)

% Does NOT check out!

%% Radiation Damping
F_radiationDamping = output.bodies(1).forceRadiationDamping;

% in out 1 DOF
F_radiationDamping5 = F_radiationDamping(:,5) + F_radiationDamping(:,1)*a.*cos(x(:,5)) - F_radiationDamping(:,3)*a.*sin(x(:,5));

% Compare to C * Xdot
C = radiationDamping(1:n,1:n);
F_C_rad = C(1,4).*a.*vel(:,5).*cos(x(:,5)) + C(5,5).*vel(:,5) + C(1,1).*vel(:,5).*a^2.*(cos(x(:,5))).^2 + C(1,5).*vel(:,5).*a.*cos(x(:,5)) + C(3,3).*vel(:,5).*a^2.*(sin(x(:,5))).^2;


figure, plot(t,F_radiationDamping5,t,F_C_rad), legend('From Output','My model')
figure, plot(t,F_radiationDamping5)
figure, plot(t,F_C_rad)

% Does NOT check out!

%% Restoring Force
F_restoring = output.bodies(1).forceRestoring;

% In our 1 DOF
F_restoring5 = F_restoring(:,5) + F_restoring(:,1)*a.*cos(x(:,5)) - F_restoring(:,3)*a.*sin(x(:,5));

% compare to K* X
K = hydrostaticstiffness(1:n,1:n);
Fk = K(5,5).*x(:,5) - K(3,3).*sin(x(:,5)).*(-a-b+a.*cos(x(:,5)));
Fk = -g*1.1*Fk;

figure, plot(t,F_restoring5,t,Fk), legend('From Output','My model')
figure, plot(t,F_restoring5)
figure, plot(t,Fk)

% Does NOT check out

%% SUM UP ALL FORCES:
SUM_Forces = output.bodies(1).forceAddedMass + output.bodies(1).forceExcitation + output.bodies(1).forceLinearDamping + output.bodies(1).forceMorisonAndViscous + output.bodies(1).forceRadiationDamping + output.bodies(1).forceRestoring;
figure, plot(t,SUM_Forces)
figure, plot(t, output.bodies(1).forceTotal)

% SUM_Forces =  output.bodies(1).forceExcitation + output.bodies(1).forceLinearDamping + output.bodies(1).forceMorisonAndViscous + output.bodies(1).forceRadiationDamping + output.bodies(1).forceRestoring;
% figure, plot(t,SUM_Forces)
% figure, plot(t, output.bodies(1).forceTotal)

