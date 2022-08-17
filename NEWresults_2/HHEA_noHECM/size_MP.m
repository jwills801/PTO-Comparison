chi = PR(2)./PR(2:end); % Frac disp of the pump/motor when connected to each rail
% The other end of the pump/motor is always tank.
V = NetRailEnergy(2:end)./PR(2:end);


% V = Q*t
    % Q = chi*w*D  => V = chi*w*D*t
% T = sum(t) => each t refers to each rail
    % T = sum(V/chi)/w/D
T = 150; % seconds
w = 2000/60; % rev/s
D_MP = sum(abs(V)./chi)/w/T; % m^3/rev


