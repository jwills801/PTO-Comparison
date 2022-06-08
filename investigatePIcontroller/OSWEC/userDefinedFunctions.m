%Example of user input MATLAB file for post processing
close all

w0 = 2*pi/waves.T;
w = body(1).hydroData.simulation_parameters.w;
w_ind = find(min(abs(w-w0)) == abs(w-w0)); % index of body BEM frequency closest to w0

% if all parameters taken in dof=5, then units will be in correct rotational form
M = body(1).momOfInertia;                                                % Mass of the system (does it need to be moment of inertia for pitch?)
mAdd = body(1).hydroData.hydro_coeffs.added_mass.all(:,:,w_ind).*simu.rho;            % Added mass in pitch
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
