%% Get decisions back out

CostMatrix = f(lam);

% The actual function
J = NaN(length(t_c),length(PRA1));
DecisionMatrix = J';
J(end,:) = zeros(size(PRA1));
for t_ind = (length(t_c)-1):-1:1
    maxSwitchingLoss = max(Eloss_A1(:, :, t_ind+1 ),[],'All','includenan') + max(Eloss_B1(:, :, t_ind+1 ),[],'All','includenan');
    [min_w_out_switches, min_ind_w_out_switches] = min(CostMatrix(t_ind+1,:) + J(t_ind+1,:));
    inds_2_check = find(   (CostMatrix(t_ind+1,:) + J(t_ind+1,:))  < (maxSwitchingLoss/t_c(2) + min_w_out_switches ));
    
    SwitchingLoss = (PRA1==PR)*Eloss_A1(:, (PRA1(min_ind_w_out_switches)==PR), t_ind+1 ) + (PRB1==PR)*Eloss_B1(:, (PRB1(min_ind_w_out_switches)==PR), t_ind+1 );

    J(t_ind,:) = min_w_out_switches + SwitchingLoss/t_c(2);
    DecisionMatrix(:,t_ind) = min_ind_w_out_switches;
    for i = inds_2_check

        % Construct Switching loss vector for this specific time index and
        % valve configuration
        % i denotes the previous valve configuration, we are deciding where
        % to go next by taking the min
        SwitchingLoss = Eloss_A1((PR==PRA1(i)), :, t_ind+1 )*(PRA1==PR)' + Eloss_B1((PR==PRB1(i)), :, t_ind+1 )*(PRB1==PR)';

        [J(t_ind,i),DecisionMatrix(i,t_ind)] = min( CostMatrix(t_ind+1,:) + J(t_ind+1,:) + SwitchingLoss/t_c(2) );
        % ELoss is divided by the time step because it is an energy, while
        % the other terms are power
    end
end

[cost, starting_ind] = min(J(1,:));



Decision_vector_c = NaN(1,length(t_c)-1);
Decision_vector_c(1) = DecisionMatrix(starting_ind,1);
for t_ind_c = 2:length(t_c)-1
    Decision_vector_c(t_ind_c) = DecisionMatrix(Decision_vector_c(t_ind_c-1),t_ind_c);
end

d_ind_c = NaN(1,length(t_c));
d_ind_c(1) = sub2ind([length(t_c),length(PRA1)],1,Decision_vector_c(1));
for t_ind_c = 2:length(t_c)
    d_ind_c(t_ind_c) = sub2ind([length(t_c),length(PRA1)],t_ind_c,Decision_vector_c(t_ind_c-1));
end

ind2stop = find(round(t_c(end),5) == round(t,5))-1;
Decision_vector = NaN(1,length(t)-1);
Decision_vector(1:spacing) = Decision_vector_c(1);
for t_ind = 1:ind2stop
    if mod(t_ind,spacing) == 1
        t_ind_c = (t_ind-1)/spacing + 1;
        Decision_vector(t_ind+spacing-1) = Decision_vector_c(t_ind_c);
    else
        Decision_vector(t_ind+spacing-1) = Decision_vector(t_ind+spacing-2);
    end
end
Decision_vector((ind2stop+1):end) = Decision_vector(ind2stop);

d_ind = NaN(1,length(t));
d_ind(1) = sub2ind([length(t),length(PRA1)],1,Decision_vector(1));
for t_ind = 2:length(t)
    d_ind(t_ind) = sub2ind([length(t),length(PRA1)],t_ind,Decision_vector(t_ind-1));
end