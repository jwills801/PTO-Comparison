clear
load hydrodata.mat
n = 6;
M = Mass(1:n,1:n) + addedmass(1:n,1:n);
C = LinearDamping(1:n,1:n) + radiationDamping(1:n,1:n);
K = hydrostaticstiffness(1:n,1:n);
Fexc = Fexc(:,1:n);


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

%% Try including all 6 DOF
LHS_6  = M*acc' + C*vel' + K*x(:,1:6)';

figure, plot(t,LHS_6,t,Fexc'), legend('Calculated','Actual'), ylabel('Force (N)'), xlabel('Time(s)')
figure, plot(t,LHS_6-Fexc'), ylabel('Force (N)'), xlabel('Time(s)')

%% Try integrating myself instead of multiplying
params = struct(); params.t = t; params.Fexc5 = Fexc5; params.M = M; params.C = C; params.K = K; params.a = a; params.b = b;
[t_,theta_] = ode23t(@(t,theta)fun(theta,t,params),0:300,[0;0]);
figure,plot(t_,theta_(:,1),t,x(:,5)), legend('My Integration','Original Data'), ylabel('Theta (rad/s)'), xlabel('Time (s)')


function thetadot = fun(theta,t_,params)
% unpack params
t= params.t; Fexc5=params.Fexc5;M= params.M;C= params.C;K= params.K;a=params.a; b=params.b;

% Derivative funtion
thetadot = [theta(2);( (-(K(5,5)*theta(1) - K(3,3)*sin(theta(1))*(-a-b+a.*cos(theta(1)))) -(C(1,4)*a*theta(2)*cos(theta(1)) + C(5,5)*theta(2) + C(1,1)*theta(2)*a^2*(cos(theta(2)))^2 + C(1,5)*theta(2)*a.*cos(theta(1)) + C(3,3)*theta(2)*a^2*(sin(theta(1)))^2 )    -(M(3,3)*a^2*sin(theta(1))*cos(theta(1)) - M(1,1)*a^2*sin(theta(1))*cos(theta(1)) - M(5,1)*a*sin(theta(1)))*theta(2)^2  + interp1(t,Fexc5,t_)) /  (M(5,1)*a*cos(theta(1)) + M(5,5) + M(1,1)*a^2*(cos(theta(1)))^2 + M(1,5)*a*cos(theta(1)) + M(3,3)*a^2*(sin(theta(1)))^2  )  ) ];
end

