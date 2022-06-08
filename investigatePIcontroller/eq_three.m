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

LHS = M*x_mddot + C*x_mdot + K*x_m

p = 1000;
for i = 1:p%length(t)
    Fexc_calc(i) = eval(subs(LHS(1),[theta2 theta2dot theta2ddot],[x(i,5) vel(i,5) acc(i,5)]));
end

%Fexc_calc = subs(LHS(5),[theta2 theta2dot theta2ddot],[x(:,5) vel(:,5) acc(:,5)]);
figure, plot(t,Fexc(:,5),t(1:p),-Fexc_calc), legend('Fexc','Fexc calc'), xlim([0 t(p)])

figure, plot(t,Fexc)
% 
% p = 1000;
% for j = 1:n
%     for i = 1:p%length(t)
%         x_calc(i) = subs(x_m(j),theta2,x(i,5));
%         vel_calc(i) = subs(x_mdot(j),[theta2 theta2dot],[x(i,5) vel(i,5)]);
%         acc_calc(i) = subs(x_mddot(j),[theta2 theta2dot theta2ddot],[x(i,5) vel(i,5) acc(i,5)]);
%     end
%     
%     figure, plot(t,x(:,j),t(1:p),x_calc), legend('x','x calc'), title([num2str(j)])
%     figure, plot(t,vel(:,j),t(1:p),vel_calc), legend('vel','vel calc'), title([num2str(j)])
%     figure, plot(t,acc(:,j),t(1:p),acc_calc), legend('acc','acc calc'), title([num2str(j)])
% end

%%
t_ind = find(t==100);
Theta = [acc(t_ind:end,5)';vel(t_ind:end,5)';x(t_ind:end,5)'];
Fexc5 = Fexc(:,5) + Fexc(:,1)*a.*cos(x(:,5)) + Fexc(:,3)*a.*sin(x(:,5));
[m c k] = Fexc5(t_ind:end)'*pinv(Theta);

%%

Kp = c
Ki = k-(2*pi/8)^2*m 
% controller.Kp = 2.02e7
% controller.Ki = 3.8302e7
lhs = m*acc(t_ind:end,5) + c*vel(t_ind:end,5) + k*x(t_ind:end,5);


figure, plot(t(t_ind:end),lhs(:),t(t_ind:end),Fexc5(t_ind:end)), legend('calc','real')
