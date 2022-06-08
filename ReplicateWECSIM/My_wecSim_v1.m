clear, close all

load Li_no_Fpto_reg.mat

t = output.wave.time;
a =5; b = 3.9;

k = 0; c = 0;

% Initialize position and velocity
x = NaN(length(t),2); Fpto = zeros(size(t));
x(1,1) = output.bodies(1).position(1,5); % Initial position in pitch
x(1,2) = output.bodies(1).velocity(1,5); % Initial velocity in pitch


% Euler's method for integration
for i = 1:length(t)-1
    
    Fpto(i+1) = k*x(i,1) + c*x(i,2);

    X = [a*sin(x(i,1));0;-a-b+a*cos(x(i,1));0;x(i,1);0]; % This is the position in 6 degrees of freedom
    XDOT = [a*x(i,2)*cos(x(i,1));0;a*x(i,2)*sin(x(i,1));0;x(i,2);0]; % This is the velocity in 6 degrees of freedom
    moment_arm = [a*cos(x(i,1)) 0 -a*sin(x(i,1)) 0 1 0]; % Converts forces in 6 DOF to the only DOF we care about: Pitch

    % Calculate xddot by first calculating all the forces
    radiation_Damping = moment_arm*body(1).hydroForce.fDamping*XDOT;
    linear_Damping = moment_arm*body(1).hydroForce.linearDamping*XDOT;
    quadratic_Damping= moment_arm*body(1).hydroForce.visDrag*(XDOT.*abs(XDOT));
    stiffness = moment_arm*body(1).hydroForce.linearHydroRestCoef*(X-[body(1).cg;0;0;0]);
    buoyancy = moment_arm(3) * (body(1).mass*simu.g - simu.rho*simu.g*body(1).dispVol);
    excitation = moment_arm*output.bodies(1).forceExcitation(i,:)';
    
    % Intertial Elements - I will seperate these based on if the term contains xddot or not
        % XDDOT = [a*diff(x,t,t)*cos(x) - a*diff(x,t)^2*sin(x);0;-a*diff(x,t,t)*sin(x) - a*diff(x,t)^2*cos(x);0;diff(x,t,t);0];
        % ^ this ^ is the XDDOT, but we want to split it up into terms that contain xddot and terms that dont
    XDDOT_no_xddot = [-a*x(i,2)^2*sin(x(i,1));0;-a*x(i,2)^2*cos(x(i,1));0;0;0];
    XDDOT_with_xddot_divided_by_xddot = [a*cos(x(i,1));0;-a*sin(x(i,1));0;1;0];
    M_matrix = diag([body(1).mass,body(1).mass,body(1).mass,[body(1).momOfInertia]]);
    Interial_no_xddot = moment_arm*body(1).hydroForce.fAddedMass*XDDOT_no_xddot + moment_arm*M_matrix*XDDOT_no_xddot;

    one_over_m = 1/ ( moment_arm*body(1).hydroForce.fAddedMass*XDDOT_with_xddot_divided_by_xddot + moment_arm*M_matrix*XDDOT_with_xddot_divided_by_xddot ) ;
    xddot = one_over_m*(excitation - radiation_Damping - quadratic_Damping - stiffness - buoyancy - linear_Damping - Interial_no_xddot - Fpto(i+1));
    

    % Use derivates to step forward in time
    x(i+1,1) = x(i,1) + x(i,2)*t(2); % Step position forward one step => this is based on the velocity
    x(i+1,2) = x(i,2) + xddot*t(2); % Step position forward one step => this is based on the acceleration
end

Power_in = Fpto'.*x(:,2);
% Position
figure, plot(t,x(:,1),t,output.bodies(1).position(:,5)), legend('My Simulation','WecSim'), ylabel('Position (m)'), xlabel('Time (s)'), title('Regular Waves')

% Velocity
figure, plot(t,x(:,2),t,output.bodies(1).velocity(:,5)), legend('My Simulation','WecSim'), ylabel('Velocity (m/s)'), xlabel('Time (s)'), title('Regular Waves')

return

%% Irregular Waves
load Li_no_Fpto_irreg.mat

t = output.wave.time;
a =5; b = 3.9;

k = 0; c = 0;

% Initialize position and velocity
x = NaN(length(t),2); Fpto = zeros(size(t));
x(1,1) = output.bodies(1).position(1,5); % Initial position in pitch
x(1,2) = output.bodies(1).velocity(1,5); % Initial velocity in pitch


% Euler's method for integration
conv_step = 1e-3; % Time step for convolution integral
for i = 1:length(t)-1
    
    Fpto(i+1) = k*x(i,1) + c*x(i,2);

    X = [a*sin(x(i,1));0;-a-b+a*cos(x(i,1));0;x(i,1);0]; % This is the position in 6 degrees of freedom
    XDOT = [a*x(i,2)*cos(x(i,1));0;a*x(i,2)*sin(x(i,1));0;x(i,2);0]; % This is the velocity in 6 degrees of freedom
    moment_arm = [a*cos(x(i,1)) 0 -a*sin(x(i,1)) 0 1 0]; % Converts forces in 6 DOF to the only DOF we care about: Pitch

    % Calculate xddot by first calculating all the forces
    linear_Damping = moment_arm*body(1).hydroForce.linearDamping*XDOT;
    quadratic_Damping= moment_arm*body(1).hydroForce.visDrag*(XDOT.*abs(XDOT));
    stiffness = moment_arm*body(1).hydroForce.linearHydroRestCoef*(X-[body(1).cg;0;0;0]);
    buoyancy = moment_arm(3) * (body(1).mass*simu.g - simu.rho*simu.g*body(1).dispVol);
    excitation = moment_arm*output.bodies(1).forceExcitation(i,:)';
    
    % RADIATION DAMPING
    % For Irregular waves, this calculation requires use of the convolution integral
    for j = 0:conv_step:t(i)
        for w = waves.w

        end
    end

    
    % Intertial Elements - I will seperate these based on if the term contains xddot or not
        % XDDOT = [a*diff(x,t,t)*cos(x) - a*diff(x,t)^2*sin(x);0;-a*diff(x,t,t)*sin(x) - a*diff(x,t)^2*cos(x);0;diff(x,t,t);0];
        % ^ this ^ is the XDDOT, but we want to split it up into terms that contain xddot and terms that dont
    XDDOT_no_xddot = [-a*x(i,2)^2*sin(x(i,1));0;-a*x(i,2)^2*cos(x(i,1));0;0;0];
    XDDOT_with_xddot_divided_by_xddot = [a*cos(x(i,1));0;-a*sin(x(i,1));0;1;0];
    M_matrix = diag([body(1).mass,body(1).mass,body(1).mass,[body(1).momOfInertia]]);
    Interial_no_xddot = moment_arm*body(1).hydroForce.fAddedMass*XDDOT_no_xddot + moment_arm*M_matrix*XDDOT_no_xddot;

    one_over_m = 1/ ( moment_arm*body(1).hydroForce.fAddedMass*XDDOT_with_xddot_divided_by_xddot + moment_arm*M_matrix*XDDOT_with_xddot_divided_by_xddot ) ;
    xddot = one_over_m*(excitation - radiation_Damping - quadratic_Damping - stiffness - buoyancy - linear_Damping - Interial_no_xddot - Fpto(i+1));
    

    % Use derivates to step forward in time
    x(i+1,1) = x(i,1) + x(i,2)*t(2); % Step position forward one step => this is based on the velocity
    x(i+1,2) = x(i,2) + xddot*t(2); % Step position forward one step => this is based on the acceleration
end

Power_in = Fpto'.*x(:,2);
% Position
figure, plot(t,x(:,1),t,output.bodies(1).position(:,5)), legend('My Simulation','WecSim'), ylabel('Position (m)'), xlabel('Time (s)'), title('Irregular Waves')

% Velocity
figure, plot(t,x(:,2),t,output.bodies(1).velocity(:,5)), legend('My Simulation','WecSim'), ylabel('Velocity (m/s)'), xlabel('Time (s)'), title('Irregular Waves')

