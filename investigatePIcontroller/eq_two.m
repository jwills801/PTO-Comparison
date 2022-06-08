clear
load hydrodata.mat
n = 6;
M = Mass(1:n,1:n) + addedmass(1:n,1:n);
C = LinearDamping(1:n,1:n) + radiationDamping(1:n,1:n);
K = hydrostaticstiffness(1:n,1:n);
Fexc = Fexc(:,1:n);


x_m = zeros(length(t),n);
x_mdot = x_m;
x_mddot = x_m;
x_m(1,:) = x(1,1:n);
dt = t(2) - t(1);

Fpto = zeros(size(Fexc));

for i = 1:length(t)-1
    % EOM
    x_mddot(i+1,:) = inv(M) * ( Fexc(i,:)' + Fpto(i,:)' - C*x_mdot(i,:)' - K*x_m(i,:)' );
    
    % Constraints
    x_mddot(i+1,2) = 0;
    x_mddot(i+1,4) = 0;
    x_mddot(i+1,6) = 0;
    
    
    %integration
    x_mdot(i+1,:) = x_mdot(i,:) + x_mddot(i+1,:)*dt;
    x_m(i+1,:) = x_m(i,:) + x_mdot(i,:)*dt;
end

figure, plot(t,x_m)
figure, plot(t,x(:,1:n))

a = 5; b = 3.9;
figure, plot(t,x(:,1),t,a*sin(x(:,5))), legend('x','x calc')
figure, plot(t,x(:,3),t,-a-b+a*cos(x(:,5))), legend('z','z calc')
