%% Used to size accumulators
% Run the onepass_general code first before using this

m3toL = 1000; %cubic meters to litres
lpm = 60000; %m^3/sec to litre/min
max_acc_vol = 0;
V_est_start_1 = 0;
V_est_start_2 = 0;
V_est_start_3 = 0;
V_est_start_1_new = 0;
V_est_start_2_new = 0;
V_est_start_3_new = 0;

motor = 1; %set motor = 1 for allowing motoring of main pump

%% Calculating the flow to be provided by the pump
%set time window 
t_wind = 5; % time in seconds for which pump will provide a constant flow rate
t_wind = t_wind/(dt); %converting time (in secs) to actual time steps
t_set = 1:t_wind:length(t); 

if length(PR) == 2
    V_act = cell(1,2);
    V_act{1} = cumsum(QR{1}(d_ind)*dt); %storing volume flow associated with rail 1
    V_act{2} = cumsum(QR{2}(d_ind)*dt); %storing volume flow associated with rail 2
    V_est_1 = cell(1,(length(t_set)-1));
    V_est_2 = cell(1,(length(t_set)-1));
    V_est_1_new = cell(1,(length(t_set)-1));
    V_est_2_new = cell(1,(length(t_set)-1));

    Q_avg = cell(1,2);
    Q_avg_new = cell(1,2);
end
    
if length(PR) == 3

    V_act = cell(1,3);
    V_act{1} = cumsum(QR{1}(d_ind)*dt); %storing volume flow associated with rail 1
    V_act{2} = cumsum(QR{2}(d_ind)*dt); %storing volume flow associated with rail 2
    V_act{3} = cumsum(QR{3}(d_ind)*dt); %storing volume flow associated with rail 3
    V_est_1 = cell(1,(length(t_set)-1)); 
    V_est_2 = cell(1,(length(t_set)-1));  
    V_est_3 = cell(1,(length(t_set)-1)); 
    V_est_1_new = cell(1,(length(t_set)-1));
    V_est_2_new = cell(1,(length(t_set)-1));
    V_est_3_new = cell(1,(length(t_set)-1));
    
    Q_avg = cell(1,3);
    Q_avg_new = cell(1,3);
end

for i = 1:(length(t_set) - 1)
    if i ~= length(t_set)-1
        start = t_set(i);
        last = t_set(i+1)-1;
    else
        start = t_set(i);
        last = t_set(i+1);
    end
        
    for j = 1:length(PR)
%         if PR(1) == 0 % use when motoring is allowed & lower is not tank
        if j == 1    
            Q_avg{j}(i) = 0;
        elseif all_options{1}(j) == '0'
            Q_avg{j}(i) = 0;
        else
            Q_avg{j}(i) = (V_act{j}(last) - V_act{j}(start))/(t_wind*dt); %calculating flow rate based on original volumes
        end
        if motor == 0
            if Q_avg{j}(i) < 0
                Q_avg{j}(i) = 0; % setting neg. flow = 0 when theres no motoring
            end
        end
    end
    if length(PR) == 2
        Q_avg_new{1}(i) = Q_avg{1}(i); %transferring orig flows to new variable
        Q_avg_new{2}(i) = Q_avg{2}(i);
    else
        Q_avg_new{1}(i) = Q_avg{1}(i); %transferring orig flows to new variable
        Q_avg_new{2}(i) = Q_avg{2}(i);
        Q_avg_new{3}(i) = Q_avg{3}(i);
    end
%% Strategy Block

%Strategy 1:
% Supplying alternately - rail starting with 0 to be decided by us
%     if mod(i,2) == 1 % odd numbered time windows will have 0 flow
%         Q_avg_new{2}(i) = 0; % this rail will be 0 at the start
%     else
%         Q_avg_new{3}(i) = 0;
%     end

% Supplying based on requirment
    if abs(Q_avg_new{2}(i)) >= abs(Q_avg_new{3}(i)) %orail with greater flow requirement is chosen for supply
        Q_avg_new{3}(i) = 0;
    else
        Q_avg_new{2}(i) = 0;
    end
    
    m = 1;
    for k = start:last
        if length(PR) == 2
            V_est_1{i}(m) = (t(k) - t(start))*Q_avg{1}(i) + V_est_start_1; %estimated volumes based on original flow rates
            V_est_2{i}(m) = (t(k) - t(start))*Q_avg{2}(i) + V_est_start_2;
        else
            V_est_1{i}(m) = (t(k) - t(start))*Q_avg{1}(i) + V_est_start_1; %estimated volumes based on original flow rates
            V_est_2{i}(m) = (t(k) - t(start))*Q_avg{2}(i) + V_est_start_2;
            V_est_3{i}(m) = (t(k) - t(start))*Q_avg{3}(i) + V_est_start_3;
        end
        m = m+1;
    end
    
    V_est_start_1 = V_est_1{i}(end);
    V_est_start_2 = V_est_2{i}(end);
    V_est_start_3 = V_est_3{i}(end);

%     if i ~= 1
%         if length(PR) == 2
%             if Q_avg{2}(i-1) ~= 0 & Q_avg_new{2}(i-1) == 0
%                 Q_avg_new{2}(i) = (V_est_2{i}(end)-V_est_2{i-1}(1))/(t_wind*dt); %compensating for the strategy
%             end
%         else
%             if Q_avg{2}(i-1) ~= 0 & Q_avg_new{2}(i-1) == 0
%                 Q_avg_new{2}(i) = (V_est_2{i}(end)-V_est_2{i-1}(1))/(t_wind*dt); %compensating for the strategy
%             end
%             if Q_avg{3}(i-1) ~= 0 & Q_avg_new{3}(i-1) == 0
%                 Q_avg_new{3}(i) = (V_est_3{i}(end)-V_est_3{i-1}(1))/(t_wind*dt);
%             end
%         end 
%     end
    
    if i ~= 1
        if abs(Q_avg_new{2}(i)) > abs(Q_avg_new{3}(i))
            temp_1 = cumsum(Q_avg{2}(1:i));
            temp_2 = cumsum(Q_avg_new{2}(1:i));
            Q_avg_new{2}(i) = temp_1(end) - temp_2(end) + Q_avg_new{2}(i);
        elseif abs(Q_avg_new{3}(i)) > abs(Q_avg_new{2}(i))
            temp_1 = cumsum(Q_avg{3}(1:i));
            temp_2 = cumsum(Q_avg_new{3}(1:i));
            Q_avg_new{3}(i) = temp_1(end) - temp_2(end) + Q_avg_new{3}(i);
        end
    end
    m = 1;
    for k = start:last
%         disp(k)
        if length(PR) == 2
            V_est_1_new{i}(m) = (t(k) - t(start))*Q_avg_new{1}(i) + V_est_start_1_new; %volumes based on compensated new flows
            V_est_2_new{i}(m) = (t(k) - t(start))*Q_avg_new{2}(i) + V_est_start_2_new;
        else
            V_est_1_new{i}(m) = (t(k) - t(start))*Q_avg_new{1}(i) + V_est_start_1_new; %volumes based on compensated new flows
            V_est_2_new{i}(m) = (t(k) - t(start))*Q_avg_new{2}(i) + V_est_start_2_new;
            V_est_3_new{i}(m) = (t(k) - t(start))*Q_avg_new{3}(i) + V_est_start_3_new;
        end
        m = m+1;
    end
    
    V_est_start_1_new = V_est_1_new{i}(end);
    V_est_start_2_new = V_est_2_new{i}(end);
    V_est_start_3_new = V_est_3_new{i}(end);
end
% return
%% Acculmulator volume:
if length(PR) == 2
    acc_max_high = abs(max(V_act{2}' - [V_est_2_new{:}]));
    acc_min_high = abs(min(V_act{2}' - [V_est_2_new{:}]));
    disp(['Maximum volume of the accumulator for rail ',num2str(length(PR))])
    acc_vol_high = (acc_max_high + acc_min_high)*m3toL
    
    acc_max_low = abs(max(V_act{1}' - [V_est_1_new{:}]));
    acc_min_low = abs(min(V_act{1}' - [V_est_1_new{:}]));
    disp(['Maximum volume of the accumulator for rail ',num2str(length(PR)-1)])
    acc_vol_low = (acc_max_low + acc_min_low)*m3toL
else
    acc_max_high = abs(max(V_act{3}' - [V_est_3_new{:}]));
    acc_min_high = abs(min(V_act{3}' - [V_est_3_new{:}]));
    disp(['Maximum volume of the accumulator for rail ',num2str(length(PR))])
    acc_vol_high = (acc_max_high + acc_min_high)*m3toL
    
    acc_max_mid = abs(max(V_act{2}' - [V_est_2_new{:}]));
    acc_min_mid = abs(min(V_act{2}' - [V_est_2_new{:}]));
    disp(['Maximum volume of the accumulator for rail ',num2str(length(PR)-1)])
    acc_vol_mid = (acc_max_mid + acc_min_mid)*m3toL
    
    acc_max_low = abs(max(V_act{1}' - [V_est_1_new{:}]));
    acc_min_low = abs(min(V_act{1}' - [V_est_1_new{:}]));
    disp(['Maximum volume of the accumulator for rail ',num2str(length(PR)-2)])
    acc_vol_low = (acc_max_low + acc_min_low)*m3toL
end
    
%maximum flow required to be provided by the pump for rail 3:
disp(['Max flow rate by main pump for rail',num2str(length(PR))])
if length(PR) == 2
    max(abs(Q_avg_new{2}))*lpm %max flow rate required for rail 2 in lpm
    disp(['Max flow rate by main pump for rail',num2str(length(PR)-1)])
    max(abs(Q_avg_new{1}))*lpm
else
    max(abs(Q_avg_new{3}))*lpm %max flow rate required for rail 3 in lpm
    disp(['Max flow rate by main pump for rail',num2str(length(PR)-1)])
    max(abs(Q_avg_new{2}))*lpm %max flow rate required for rail 2 in lpm
    disp(['Max flow rate by main pump for rail',num2str(length(PR)-2)])
    max(abs(Q_avg_new{1}))*lpm %max flow rate required for rail 2 in lpm
end


%% Plotting Actual Flows vs the Pump Flows for each rail

figure(20)
if length(PR) == 2
    plot(t,V_act{1},'b', t,V_act{2},'b')
    hold on
    for ii = 1:(length(t_set)-1)
        if ii ~= length(t_set) - 1
            plot((t_set(ii):t_set(ii+1))*dt, V_est_1_new{ii},'r',(t_set(ii):t_set(ii+1))*dt, V_est_2_new{ii},'r')
            hold on
        else
            plot((t_set(ii):t_set(ii+1)+1)*dt, V_est_1_new{ii},'r',(t_set(ii):t_set(ii+1)+1)*dt, V_est_2_new{ii},'r')
            hold on
        end 
    end
elseif length(PR) == 3
    plot(t,V_act{1},'b',t,V_act{2},'b',t,V_act{3},'b')
    hold on
    for ii = 1:(length(t_set)-1)
        if ii ~= length(t_set) - 1
            plot((t_set(ii):t_set(ii+1)-1)*dt, V_est_1_new{ii}, 'r', (t_set(ii):t_set(ii+1)-1)*dt, V_est_2_new{ii}, 'r',...
                (t_set(ii):t_set(ii+1)-1)*dt, V_est_3_new{ii}, 'r')
            hold on
        else
            plot((t_set(ii):t_set(ii+1))*dt, V_est_1_new{ii}, 'r', (t_set(ii):t_set(ii+1))*dt, V_est_2_new{ii}, 'r',...
                (t_set(ii):t_set(ii+1))*dt, V_est_3_new{ii}, 'r')
            hold on
        end
    end    
end
xlabel('Time in seconds')
ylabel('Volume in m^3')
title('Actual flow requirements vs average pump supply')
% legend('Required Flow','Provided Flow')
    
        
        
    
    
    