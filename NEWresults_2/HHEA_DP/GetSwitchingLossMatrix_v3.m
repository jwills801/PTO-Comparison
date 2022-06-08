% In this version, I will make a tabulated grid of points and just
% interpolate to find the energy loss

% Set the volume and velocity*Area parameters
velA_min = min(V1)*ACap1;
velA_max = max(V1)*ACap1;
velA_n = 10;
velA_vals = linspace(velA_min,velA_max,velA_n);

vol_min = min([min(Vol_A1),min(Vol_B1)]);
vol_max = max([max(Vol_A1),max(Vol_B1)]);
vol_n = 1;
vol_vals = linspace(vol_min,vol_max,vol_n);

[velA_matrix,vol_matrix] = ndgrid(velA_vals,vol_vals);
velA_vec = velA_matrix(:); vol_vec = vol_matrix(:);
Eloss = NaN(length(PR),length(PR),length(velA_vec)); delaychosen = NaN(length(PR),length(PR),length(velA_vec));

PA = NaN(size(tspan')); Q_C = NaN(length(tspan)-1,1); Q_O = NaN(length(tspan)-1,1);
Eloss_delay = NaN(size(delayvals));
for i = 1:length(velA_vec)
    for j = 1:length(PR)
        for k = 1:length(PR)
            if j == k
                Eloss(j,k,i) = (abs(velA_vec(i)))^3/k_/k_*Tf;
                delaychosen(j,k,i) = NaN;
            else
                Eloss_delay(:) = NaN;
                for ii = 1:length(delayvals)
                    delay = round(delayvals(ii));
                    
                    
                    xon = step(gs,tspan);
                    xoff = (1-xon);
                    xon = [zeros(delay,1);xon(1:end-delay)];
                    
                    % The is no oscillation in the vavle. once it opens, it stays open
                    % Once it closes it closes
                    xoff = max(0,xoff);
                    xon = min(1,xon);
                    % Also, once it has reaches open or close, dont let it go back at all
                    for jj = 5:length(tspan)
                        if xoff(jj-1) == 0
                            xoff(jj) = 0;
                        elseif xon(jj-1) == 1
                            xon(jj) = 1;
                        end
                    end
                    
                    PA(1) = PR(j);
                    for kk=1:length(tspan)-1
                        Q_C(kk) = k_*xoff(kk)*sign(PR(j)-PA(kk)).*sqrt(abs(PR(j)-PA(kk))); % Oriface equation for closing valve
                        Q_O(kk) = k_*xon(kk)*sign(PR(k)-PA(kk)).*sqrt(abs(PR(k)-PA(kk))); % Oriface equation for opening valve
                        PA(kk+1) = PA(kk) + beta/vol_vec(i)*(Q_C(kk)+Q_O(kk)-velA_vec(i))*dt; % Compressibility -> expression for dPdt
                    end
                    Loss_h = Q_C.*(PR(j)-PA(1:length(tspan)-1));
                    Loss_l = Q_O.*(PR(k)-PA(1:length(tspan)-1));
                    Eloss_delay(ii) =  sum(Loss_h+Loss_l)*dt;
                end
                [Eloss(j,k,i), delaychosen(j,k,i)] = min(Eloss_delay);
            end
            if t(2) > Tf
                Eloss(j,k,i) = Eloss(j,k,i)+(abs(velA_vec(i)))^3/k_/k_*(t(2)-Tf);
            end
        end
    end
end

%% Now map this to the drive cycle
Eloss_A1 = NaN(length(PR),length(PR),length(t)); Eloss_B1 = NaN(length(PR),length(PR),length(t)); Eloss_A2 = NaN(length(PR),length(PR),length(t)); Eloss_B2 = NaN(length(PR),length(PR),length(t));
for i = 1:length(t)
    
[~,vol_ind_A1] = min(abs(vol_vals-Vol_A1));
[~,velA_ind_A1] = min(abs(velA_vals-V1(i)*ACap1));
Eloss_A1(:,:,i) = Eloss(:,:,sub2ind(size(velA_matrix),velA_ind_A1,vol_ind_A1));

[~,vol_ind_B1] = min(abs(vol_vals-Vol_B1));
[~,velA_ind_B1] = min(abs(velA_vals+V1(i)*ARod1));
Eloss_B1(:,:,i) = Eloss(:,:,sub2ind(size(velA_matrix),velA_ind_B1,vol_ind_B1));
end


% delaychosen(:,:,sub2ind(size(velA_matrix),velA_ind_A1,vol_ind_A1)); %
% Which delay is used may be useful later, but right now it is just tedious

%% Putting them all together into one matrix
Eloss_SwitchingValve = NaN(length(PRA1),length(PRA1),length(t));
for i = 1:length(PRA1)
    for j = 1:length(PRA1)
        ind_A1_old = find(PRA1(i) == PR); ind_B1_old = find(PRB1(i) == PR);
        ind_A1_new = find(PRA1(j) == PR); ind_B1_new = find(PRB1(j) == PR);
        
        Eloss_SwitchingValve(i,j,:) = Eloss_A1(ind_A1_old,ind_A1_new,:)  + Eloss_B1(ind_B1_old,ind_B1_new,:);
    end
end

