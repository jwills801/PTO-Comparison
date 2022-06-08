%% Non-Compressible Fluid Hydraulic PTO-Sim

ptosim = ptoSimClass('Non-Compressible Fluid Hydraulic');

% Rotary to Linear Adjustable Rod
% same sizing as PTO-Sim OSWEC example
ptosim.motionMechanism.crank = 3;
ptosim.motionMechanism.offset = 1.3;
ptosim.motionMechanism.rodInit = 5;

%% Controller Parameters
% This portion of the input file defines the PI gains for the controller.
% Based on SAND2018-10945. See section 1.4.1
controller = struct();
controller.name = 'PI';

% sum refers to using PI(u_actual-u_opt) or just PI(u_actual).
%  u_optimal = F_exc/(2*Ri)
% TODO Which to do? --> pretty sure its supposed to be 0, just take actual velocity
controller.pi_sum = 0;

w0 = 2*pi/waves.T;
w = body(1).hydroData.simulation_parameters.w;
w_ind = find(min(abs(w-w0)) == abs(w-w0)); % index of body BEM frequency closest to w0

% if all parameters taken in dof=5, then units will be in correct rotational form
M = body(1).momOfInertia(2);                                                % Mass of the system (does it need to be moment of inertia for pitch?)
mAdd = body(1).hydroData.hydro_coeffs.added_mass.all(5,5,w_ind).*simu.rho;            % Added mass in pitch
R = body(1).hydroData.hydro_coeffs.radiation_damping.all(5,5,w_ind).*simu.rho*2*pi/waves.T;        % Radiation damping in pitch
Hs_s = body(1).hydroData.hydro_coeffs.linear_restoring_stiffness(5,5)*simu.rho*simu.g;      % Linear, hydrostatic restoring coefficient in pitch

% From paper: Bv =  'linear damping term describing the viscous
% effects of the fluid (and/or any other linear friction terms)'
%     -> include all damping effects except radiation
% Bv = body linearDamping + PTO damping + mooring damping + body viscousDamping + ...
% Bv + body(1).linearDamping + pto(1).c + mooring(1).c + body(1).hydroForce.visDrag*|Velocity(5)|
%     body(1).hydroForce.visDrag = diag(0.5*rho.*body(1).viscDrag.cd.*body(1).viscDrag.characteristicArea);
viscLDamp = body(1).linearDamping(5,5); % body linear damping
ptoDamp = pto(1).c; % PTO linear damping (rotational for rotational pto by default)
cd = body(1).viscDrag.cd(5); % body quadratic viscous drag coefficient
ca = body(1).viscDrag.characteristicArea(5); % body quadratic viscous drag area
viscQDamp = diag(0.5*simu.rho*cd*ca)*w0; % body quadratic viscous drag, F_visc / vel = c_visc
Bv = viscLDamp + ptoDamp + viscQDamp;                                      % Linear damping coefficient in pitch

controller.Kp = Bv + R;                                           % Proportional gain (damping)
controller.Ki = w0^2*(M+mAdd) + Hs_s;                            % Integral gain (mass, stiffness)

%clear w0 w w_ind M mAdd Bv R Hs_s hydro viscLDamp ptoDamp viscQDamp

controller.Kp = 2.02e7;
controller.Ki = 3.8302e7;


controller.Kp = 0;
controller.Ki = 0;

% Old values:
% regular wave: Kp 1e7, Ki -9.4662e4
% irregular wave: Kp 1e7, Ki -1.5445e6