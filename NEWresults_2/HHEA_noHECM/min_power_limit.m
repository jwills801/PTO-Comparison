function [Max_batt_pow] = min_power_limit(battery_power,P1_Cyl,t_c,spacing)
% This code is called in Make_losses_DP
% It finds the smallest possible torque limit for each actuator.
% This is important because as the time between switches increases, the torque requirment on the HECM may become quite large.

% inputs:
% T1_Act: a matrix with the torque on the HECM shaft if a certain rail configuration is chosen at a certain time.
% rows denote time, columns denote rail configurations
% P1_Cyl: a matrix with the pressure in the rod side (the outlet of the HECM) if a certain rail configuration is chosen at a certain time.
% rows denote time, columns denote rail configurations
% t_c is a uniform vector of time points at which decisions can be made.
% spacing is the number of fine time steps (t) inside of one decision time step (t_c).

min_Pow = 0;
important_time_ind_c = 0;

% Cavitation constraint
battery_power(P1_Cyl < -1e5) = inf;
for t_ind_c = 1:length(t_c)-1
    t_ind = (t_ind_c-1)*spacing+1;
    
    if spacing == 1
        temp = min(abs(battery_power(t_ind:(t_ind+spacing-1),:))); % min power at this time interval
    else
        temp = min(max(abs(battery_power(t_ind:(t_ind+spacing-1),:)))); % min power at this time interval
    end
    if temp == inf
        error(['Homemade Error: Infesibility occurs at time ' num2str(t_c(t_ind_c)) ' seconds'])
    end
    [min_Pow, temp_ind] = max([min_Pow,temp]);
    if temp_ind == 2
        important_time_ind_c = t_ind_c;
    end
end
disp(['Power constraint: ' num2str(round(1.01*min_Pow/1000,3)) ' kW. Set at time interval begining at ' num2str(t_c(important_time_ind_c)) ' seconds'])
Max_batt_pow = 1.01*min_Pow; % Give a one percent breathing room
end