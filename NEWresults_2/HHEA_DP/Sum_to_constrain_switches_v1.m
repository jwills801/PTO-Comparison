%% Constrain Time between switches
tic
%Preallocate memory:
HECMLosses_c = NaN(length(t_c),length(PRA1)); battery_power_c = NaN(length(t_c),length(PRA1)); 
for i = 1:length(PR)
    QR_c{i} = NaN(length(t_c),length(PRA1));
end

for t_ind_c = 1:length(t_c)
    if length(t_c) == length(t)
        HECMLosses_c = HECMLosses;
        battery_power_c = battery_power;
    else
        if t_ind_c == length(t_c)
            t_ind_low = find(round(t_c(t_ind_c),3) == round(t,3));
            if isempty(t_ind_low)
                disp('Something went wrong while summing')
            end

            HECMLosses_c(t_ind_c,:) = sum(HECMLosses(t_ind_low:end,:),1) ;
            battery_power_c(t_ind_c,:) = sum(battery_power(t_ind_low:end,:),1) ;
            for i = 1:length(PR)
                QR_c{i}(t_ind_c,:) = sum(QR{i}(t_ind_low:end,:),1);
            end
        else
            t_ind_low = (t_ind_c-1)*spacing+1;
            t_ind_high = (t_ind_c)*spacing+1;

            if isempty(t_ind_low) | isempty(t_ind_high)
                disp('Something went wrong while summing')
            end

            HECMLosses_c(t_ind_c,:) = sum(HECMLosses(t_ind_low:t_ind_high-1,:),1) ;
            battery_power_c(t_ind_c,:) = sum(battery_power(t_ind_low:t_ind_high-1,:),1) ;
            for i = 1:length(PR)
                QR_c{i}(t_ind_c,:) = sum(QR{i}(t_ind_low:t_ind_high-1,:),1);
            end
        end
    end
end
disp(['Summing to constrain swtiches took ' num2str(toc) ' seconds'])
% Check fesibility
feasibility = isfinite(  max(min(HECMLosses_c + battery_power_c,[],2),[],'includenan') )

if feasibility == 0
    [~,problem_ind_c] = max(min(HECMLosses_c + battery_power_c,[],2));

    figure, plot(t_c, HECMLosses_c), xlabel('Time (s)'), ylabel('HECM losses (W)')

    error(['Homemade Error: infeasiblity at time ' num2str(t_c(problem_ind_c)) ' seconds'])
end
