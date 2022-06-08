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


t_ind = find(t==100.44);
Theta = [acc(t_ind:end,5)';vel(t_ind:end,5)';vel(t_ind:end,5)'.^2;x(t_ind:end,5)'];
Fexc5 = Fexc(:,5) + Fexc(:,1)*a.*cos(x(:,5)) - Fexc(:,3)*a.*sin(x(:,5));
params = Fexc5(t_ind:end)'*pinv(Theta);

m = params(1)
c = params(2)
k = params(3)



Kp = c
Ki = k-(2*pi/8)^2*m 
% controller.Kp = 2.02e7
% controller.Ki = 3.8302e7
lhs = params*Theta;


figure, plot(t(t_ind:end),lhs(:),t(t_ind:end),Fexc5(t_ind:end)), legend('calc','real')


%%

lhs2 = M(5,1)*acc(:,1) + M(5,5)*acc(:,5) + C(5,1)*vel(:,1) + C(5,5)*vel(:,5) + K(5,5)*x(:,5);
figure, plot(t,lhs2,t,Fexc5), legend('calc','real')