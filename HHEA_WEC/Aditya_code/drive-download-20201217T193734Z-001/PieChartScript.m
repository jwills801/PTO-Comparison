%% Pie Chart For HHEA Script
close all

    %% Load Sensing general - SELECT FOR COMPARISON WITH BASIC HHEA
    OutL_JCB05_ND = [485.6311  264.8924  801.4999]; % LS - ND
    InsL_JCB05_ND = 1e3*    [0.2081    1.3440];
    
    OutL_JCB05_GR = 1e3*    [0.3790    0.2384    1.1870]; % LS - GR 
    InsL_JCB05_GR = 1e3*    [0.1418    1.6624];
    
    OutL_JCB05_TR = [436.2816  243.0339  668.0590]; % LS - TR
    InsL_JCB05_TR = 1e3*    [0.0840    1.2634];
    
    %% Load Sensing with battery - SELECT FOR COMPARISON WITH FE-HHEA
%     OutL_JCB05_ND = [485.6311  289.4162  801.4999]; % LS - ND
%     InsL_JCB05_ND = 1e3*    [0.2081    1.3685];
%     
%     OutL_JCB05_GR = 1e3*    [0.3790    0.2722    1.1870]; % LS - GR 
%     InsL_JCB05_GR = 1e3*    [0.1418    1.6963];
%     
%     OutL_JCB05_TR = [436.2816  271.6854  668.0590]; % LS - TR
%     InsL_JCB05_TR = 1e3*    [0.0840    1.2921];
    
    %% Values for 2 CPR Case - FE-HHEA (Optimal Solution = NP)
    OutH2_JCB05_ND = [485.6311   32.1886  172.1647]; % HHEA - ND 
    InsH2_JCB05_ND = [208.0657  481.9188];
        
    OutH2_JCB05_GR = [379.4082  115.7419  318.4305]; % HHEA - GR 
    InsH2_JCB05_GR = [142.7560  670.8246];
        
    OutH2_JCB05_TR = [436.6190   27.2265  117.9133]; % STEAM - GR
    InsH2_JCB05_TR = [84.2683  497.4905];
    
    %% Values for 2 CPR Case - FE-HHEA (Optimal Solution = C0)
%     OutH2_NP_JCB05_ND = [485.6311    0.0073  165.2109]; % HHEA - ND
%     InsH2_NP_JCB05_ND = [208.0657  442.7837]; % for [70 70 70 43]

    OutH2_NP_JCB05_ND = [485.6311    0.0011  185.2925]; % HHEA - ND  
    InsH2_NP_JCB05_ND = [208.0657  462.8591]; % for [24 78 39 43]

%     OutH2_NP_JCB05_GR = [379.4082    0.0018  370.6666]; % HHEA - GR 
%     InsH2_NP_JCB05_GR = [142.7560  607.3206];% for [70 70 70 43]

    OutH2_NP_JCB05_GR = [379.4082    0.0017  355.7732]; % HHEA - GR 
    InsH2_NP_JCB05_GR = [142.7560  592.4270]; % for [24 78 39 43]
        
%     OutH2_NP_JCB05_TR = [436.6190    0.0046  115.6430]; % STEAM - GR
%     InsH2_NP_JCB05_TR = [84.2683  467.9983];% for [70 70 70 43]

    OutH2_NP_JCB05_TR = [436.6190    0.0011  132.8488]; % STEAM - GR
    InsH2_NP_JCB05_TR = [84.2683  485.2006];% for [24 78 39 43]

    %% Values for 3 CPR Case - FE-HHEA (Optimal Solution = NPP)
    OutH3_JCB05_ND = [485.6311   40.9302  137.5935]; % HHEA - ND 
    InsH3_JCB05_ND = [208.0657  456.0892];
        
    OutH3_JCB05_GR = [379.4082   58.2989  268.9312]; % HHEA - GR 
    InsH3_JCB05_GR = [142.7560  563.8823];
    
    OutH3_JCB05_GR_Base = [379.4082 58.69 233.45]; % Basic_HHEA - GR - 'C0P'
    InsH3_JCB05_GR_Base = [142.76 528.80];
        
    OutH3_JCB05_TR = [436.6190   50.9460   97.2466]; % STEAM - GR
    InsH3_JCB05_TR = [84.2683  500.5434];
    
    %% Values for 3 CPR Case - FE-HHEA (Optimal Solution = C00)
%     OutH3_NP_JCB05_ND = [485.6311    0.0168  151.9821]; % HHEA - ND - Constrained sizes
%     InsH3_NP_JCB05_ND = [208.0657  429.5643]; 
    
    OutH3_NP_JCB05_ND = []; % HHEA - ND - Unconstrained sizes
    InsH3_NP_JCB05_ND = []; 
        
%     OutH3_NP_JCB05_GR = [379.4082    0.0197  258.9473]; % HHEA - GR - Constrained sizes
%     InsH3_NP_JCB05_GR = [142.7560  495.6192];
    
    OutH3_NP_JCB05_GR = [379.4082    0.0053  229.1277]; % HHEA - GR - Unconstrained sizes
    InsH3_NP_JCB05_GR = [142.7560  465.7852];
        
%     OutH3_NP_JCB05_TR = [436.6190    0.0013  121.0180]; % STEAM - GR - Constrained sizes
%     InsH3_NP_JCB05_TR = [84.2683  473.3700];

    OutH3_NP_JCB05_TR = []; % STEAM - GR - Unconstrained sizes
    InsH3_NP_JCB05_TR = [];
    

% pie_array = [PositiveWork MainPumpLosses  ControlLosses]
%cake_array = [NegativeWork Input_Energy]
%% Pie Charts fo the 2 CPR Case: Optimal Solution = NP
%     pie_array = OutH2_JCB05_ND/sum(InsL_JCB05_ND) % Sim 1 HHEA
%     cake_array = InsH2_JCB05_ND/sum(InsL_JCB05_ND)

%     pie_array = OutH2_JCB05_GR/sum(InsL_JCB05_GR) % Sim 2 HHEA
%     cake_array = InsH2_JCB05_GR/sum(InsL_JCB05_GR)

%     pie_array = OutH2_JCB05_TR/sum(InsL_JCB05_TR) % Sim 3 HHEA
%     cake_array = InsH2_JCB05_TR/sum(InsL_JCB05_TR)
    
%% Pie Charts for the 2 CPR case: Optimal Solution = C0
%     pie_array = OutH2_NP_JCB05_ND/sum(InsL_JCB05_ND) % Sim 1 HHEA
%     cake_array = InsH2_NP_JCB05_ND/sum(InsL_JCB05_ND)

%     pie_array = OutH2_NP_JCB05_GR/sum(InsL_JCB05_GR) % Sim 2 HHEA
%     cake_array = InsH2_NP_JCB05_GR/sum(InsL_JCB05_GR)
% % 
%     pie_array = OutH2_NP_JCB05_TR/sum(InsL_JCB05_TR) % Sim 3 HHEA
%     cake_array = InsH2_NP_JCB05_TR/sum(InsL_JCB05_TR)

%% Pie charts for the 3 CPR case: Optimal Solution = NPP
%     pie_array = OutH3_JCB05_ND/sum(InsL_JCB05_ND) % Sim 1 HHEA
%     cake_array = InsH3_JCB05_ND/sum(InsL_JCB05_ND)

%     pie_array = OutH3_JCB05_GR/sum(InsL_JCB05_GR) % Sim 2 HHEA
%     cake_array = InsH3_JCB05_GR/sum(InsL_JCB05_GR)

%     pie_array = OutH3_JCB05_GR_Base/sum(InsL_JCB05_GR) % Sim 2 HHEA
%     cake_array = InsH3_JCB05_GR_Base/sum(InsL_JCB05_GR)

%     pie_array = OutH3_JCB05_TR/sum(InsL_JCB05_TR) % Sim 3 HHEA
%     cake_array = InsH3_JCB05_TR/sum(InsL_JCB05_TR)

%% Pie Charts for the 3 CPR case: Optimal Solution = C00
%     pie_array = OutH3_NP_JCB05_ND/sum(InsL_JCB05_ND) % Sim 1 HHEA
%     cake_array = InsH3_NP_JCB05_ND/sum(InsL_JCB05_ND)

    pie_array = OutH3_NP_JCB05_GR/sum(InsL_JCB05_GR) % Sim 2 HHEA
    cake_array = InsH3_NP_JCB05_GR/sum(InsL_JCB05_GR)

%     pie_array = OutH3_NP_JCB05_TR/sum(InsL_JCB05_TR) % Sim 3 HHEA
%     cake_array = InsH3_NP_JCB05_TR/sum(InsL_JCB05_TR)
%% Pie Charts for the Load Sensing Pump 
% pie_array = OutL_JCB05_ND/sum(InsL_JCB05_ND) % Sim 1 LS
% cake_array = InsL_JCB05_ND/sum(InsL_JCB05_ND)
% 
% pie_array = OutL_JCB05_GR/sum(InsL_JCB05_GR) % Sim 2 LS
% cake_array = InsL_JCB05_GR/sum(InsL_JCB05_GR)
% % 
% pie_array = OutL_JCB05_TR/sum(InsL_JCB05_TR) % Sim 3 LS
% cake_array = InsL_JCB05_TR/sum(InsL_JCB05_TR)



%% Plot the pie chart
figure(13)
p = pie(pie_array)
pText = findobj(p,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Positive Work: ','Main Pump Losses: ','HECM Losses: '}; 
combinedtxt = strcat(txt,percentValues'); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);

%%
delete(findobj(p,'Type','text'))%option 2

