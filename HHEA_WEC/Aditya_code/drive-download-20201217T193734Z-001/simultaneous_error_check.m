%% initialize
varsbefore = who;
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
t_sample = 100; %sampling time = ((1/100)*t_sample) seconds.

t_set = t(1:t_sample:end); %sampling time reduction

motor = 1; %set motor = 1 for allowing motoring of main pump

% Q_P = ones(1,3)*20/lpm;
% Q_set = 1:20; %flow rate in lpm
Q_set = 8; %flow rate in lpm (uncomment and change to optimal Q value found by using the above set)
error_tol = 1e-4;
%% Calculating the flow to be provided by the pump

for jj = 1:length(Q_set)
    Q_lim = Q_set(jj);
    
    if length(PR) == 3

        V_act = cell(1,3);
        error = cell(1,3);
        V_act{1} = cumsum(QR{1}(d_ind)*dt); %storing volume flow associated with rail 1
        V_act{2} = cumsum(QR{2}(d_ind)*dt); %storing volume flow associated with rail 2
        V_act{3} = cumsum(QR{3}(d_ind)*dt); %storing volume flow associated with rail 3
        V_act{1} = V_act{1}(1:t_sample:end);
        V_act{2} = V_act{2}(1:t_sample:end);
        V_act{3} = V_act{3}(1:t_sample:end);
        

    end

    error_0 = 0;
%     t_0 = t(1);
    t_0 = t_set(1) %updated sampling rate
    
    V_init_2 = V_act{2}(1);
    V_init_3 = V_act{3}(1);
    V_test_2(1) = V_init_2;
    V_test_3(1) = V_init_3;

    error{2}(1) = V_act{2}(1) - V_test_2(1);
    error{3}(1) = V_act{3}(1) - V_test_3(1);

    pump = 2;
    Q_P(2) = Q_lim/lpm;
    Q_P(3) = 0; %starting with no flow to rail 3

    % pump = 3;
    % Q_P(3) = Q_lim/lpm;
    % Q_P(2) = 0; %starting with no flow to rail 2 


%     for i = 2:length(t)
    for i = 2:length(t_set) %updated sampling rate
        
        V_test_1(i) = 0;

        if all_options{1}(2) == '0'
            V_test_2(i) = 0;
        else
%             V_test_2(i) = Q_P(2)*(t(i)-t_0) + V_init_2;
            V_test_2(i) = Q_P(2)*(t_set(i)-t_0) + V_init_2; %updated sampling rate
            error{2}(i) = V_act{2}(i) - V_test_2(i);
        end

        if all_options{1}(3) == '0'
            V_test_3(i) = 0;
        else
%              V_test_3(i) = Q_P(3)*(t(i)-t_0) + V_init_3;
             V_test_3(i) = Q_P(3)*(t_set(i)-t_0) + V_init_3; %updated sampling rate
             error{3}(i) = V_act{3}(i) - V_test_3(i);
        end
        
        if error{2}(i) >= error{3}(i)
            Q_P(2) = Q_lim/lpm;
            Q_P(3) = 0;
        else
            Q_P(2) = 0;
            Q_P(3) = Q_lim/lpm;
        end

        
        

    %     V_test_2(i) = Q_P(2)*(t(i)-t_0) + V_init_2;
    %     error{2}(i) = V_act{2}(i) - V_test_2(i);
    %     
    %     V_test_3(i) = Q_P(3)*(t(i)-t_0) + V_init_3;
    %     error{3}(i) = V_act{3}(i) - V_test_3(i);

        Q_P_2_set(i) = Q_P(2);
        Q_P_3_set(i) = Q_P(3);

        V_init_2 = V_test_2(i);
        V_init_3 = V_test_3(i);

        t_0 = t_set(i);
    end

    acc_size_2 = (abs(max(V_act{2} - V_test_2')) + abs(min(V_act{2} - V_test_2')))*m3toL;

    acc_size_3 = (abs(max(V_act{3} - V_test_3')) + abs(min(V_act{3} - V_test_3')))*m3toL;
    
    acc_size_2_set(jj) = acc_size_2;
    acc_size_3_set(jj) = acc_size_3;
end
idx = find((acc_size_2_set + acc_size_3_set) == min(acc_size_2_set + acc_size_3_set))

acc_size_2 = acc_size_2_set(idx)
acc_size_3 = acc_size_3_set(idx)
Q_lim = Q_set(idx)

%% plot for that flow rate

figure(24)
% plot(t,V_test_1,'r',t,V_test_2,'r',t,V_test_3,'r')
plot(t_set,V_test_1,'r',t_set,V_test_2,'r',t_set,V_test_3,'r') %updated sampling rate
hold on
% plot(t,V_act{1},'b',t,V_act{2},'b',t,V_act{3},'b')
plot(t_set,V_act{1},'b',t_set,V_act{2},'b',t_set,V_act{3},'b') %updated sampling rate
xlabel('Time in seconds')
ylabel('Cumulative Vol. (m^3)')
% legend('Actual Rail 2 flow', 'Actual Rail 3 flow', 'Desired Rail 2 flow', 'Desired Rail 3 flow')
title(['Flow comparison for Q =',num2str(Q_lim),' and error tolerance = ',num2str(error_tol),' m^3'])

return
%% post code run
varsafter = [];
varsnew = [];
varsafter = who;
varsnew = setdiff(varsafter,varsbefore);
clear(varsnew{:})

