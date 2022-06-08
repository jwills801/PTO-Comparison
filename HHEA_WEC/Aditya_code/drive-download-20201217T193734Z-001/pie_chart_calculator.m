    %% Load Sensing general - SELECT FOR COMPARISON WITH BASIC HHEA
%     OutL_JCB05_ND = [485.6311  264.8924  801.4999]; % LS - ND
%     InsL_JCB05_ND = 1e3*    [0.2081    1.3440];
%     
%     OutL_JCB05_GR = 1e3*    [0.3790    0.2384    1.1870]; % LS - GR 
%     InsL_JCB05_GR = 1e3*    [0.1418    1.6624];
%     
%     OutL_JCB05_TR = [436.2816  243.0339  668.0590]; % LS - TR
%     InsL_JCB05_TR = 1e3*    [0.0840    1.2634];
%     
    %% Load Sensing with Battery - SELECT FOR FE-HHEA
    OutL_JCB05_ND = [485.6311  289.4162  801.4999]; % LS - ND
    InsL_JCB05_ND = 1e3*    [0.2081    1.3685];
    
    OutL_JCB05_GR = 1e3*    [0.3790    0.2722    1.1870]; % LS - GR 
    InsL_JCB05_GR = 1e3*    [0.1418    1.6963];
    
    OutL_JCB05_TR = [436.2816  271.6854  668.0590]; % LS - TR
    InsL_JCB05_TR = 1e3*    [0.0840    1.2921];
    
    %% Values for 2 CPR Case (Optimal Solution = NP)
    OutH2_JCB05_ND = [485.6311   32.1886  172.1647]; % HHEA - ND 
    InsH2_JCB05_ND = [208.0657  481.9188];
        
    OutH2_JCB05_GR = [379.4082  115.7419  318.4305]; % HHEA - GR 
    InsH2_JCB05_GR = [142.7560  670.8246];
        
    OutH2_JCB05_TR = [436.6190   27.2265  117.9133]; % STEAM - GR
    InsH2_JCB05_TR = [84.2683  497.4905];
    
    %% Values for 2 CPR Case (Optimal Solution = C0)
%     OutH2_NP_JCB05_ND = [485.6311    0.0073  165.2109]; % HHEA - ND 
%     InsH2_NP_JCB05_ND = [208.0657  442.7837];
    
    OutH2_NP_JCB05_ND = [485.6311    0.0011  185.2925]; % HHEA - ND  
    InsH2_NP_JCB05_ND = [208.0657  462.8591]; % for [24 78 39 43
         
%     OutH2_NP_JCB05_GR = [379.4082    0.0018  370.6666]; % HHEA - GR 
%     InsH2_NP_JCB05_GR = [142.7560  607.3206];

    OutH2_NP_JCB05_GR = [379.4082    0.0017  355.7732]; % HHEA - GR 
    InsH2_NP_JCB05_GR = [142.7560  592.4270]; % for [24 78 39 43]
        
%     OutH2_NP_JCB05_TR = [436.6190    0.0046  115.6430]; % STEAM - GR
%     InsH2_NP_JCB05_TR = [84.2683  467.9983];

    OutH2_NP_JCB05_TR = [436.6190    0.0011  132.8488]; % STEAM - GR
    InsH2_NP_JCB05_TR = [84.2683  485.2006];% for [24 78 39 43]
    %% Values for 3 CPR Case (Optimal Solution = NPP)
    OutH3_JCB05_ND = [485.6311   40.9302  137.5935]; % HHEA - ND 
    InsH3_JCB05_ND = [208.0657  456.0892];
        
%     OutH3_JCB05_GR = [379.4082   58.2989  268.9312]; % HHEA - GR 
%     InsH3_JCB05_GR = [142.7560  563.8823];
    
    OutH3_JCB05_GR = [379.4082 58.69 233.45]; % Basic_HHEA - GR 
    InsH3_JCB05_GR = [142.76 528.80];
        
    OutH3_JCB05_TR = [436.6190   50.9460   97.2466]; % STEAM - GR
    InsH3_JCB05_TR = [84.2683  500.5434];
    
    %% Values for 3 CPR Case (Optimal Solution = C00)
%     OutH3_NP_JCB05_ND = [485.6311    0.0168  151.9821]; % HHEA - ND - constrained Sizes
%     InsH3_NP_JCB05_ND = [208.0657  429.5643];
    
    OutH3_NP_JCB05_ND = []; % HHEA - ND - Unconstrained Sizes
    InsH3_NP_JCB05_ND = [];
        
%     OutH3_NP_JCB05_GR = [379.4082    0.0197  258.9473]; % HHEA - GR - constrained Sizes
%     InsH3_NP_JCB05_GR = [142.7560  495.6192];
    
    OutH3_NP_JCB05_GR = [379.4082    0.0053  229.1277]; % HHEA - GR - Unconstrained Sizes
    InsH3_NP_JCB05_GR = [142.7560  465.7852];
        
%     OutH3_NP_JCB05_TR = [436.6190    0.0013  121.0180]; % STEAM - GR - constrained Sizes
%     InsH3_NP_JCB05_TR = [84.2683  473.3700];

    OutH3_NP_JCB05_TR = []; % STEAM - GR - Unconstrained Sizes
    InsH3_NP_JCB05_TR = [];

    return
    %% For 2 CPR - C0 - Ninety Degree:
    Energy_Saving_ND = sum(InsL_JCB05_ND) - sum(InsH2_NP_JCB05_ND)
    Energy_Saving_ND_Percent = 100*Energy_Saving_ND/sum(InsL_JCB05_ND)
    
    Regen_Energy_ND = InsH2_NP_JCB05_ND(1)
    Regen_Energy_ND_Percent = 100*InsH2_NP_JCB05_ND(1)/sum(InsL_JCB05_ND)
    Regen_Energy_ND_Pi_Angle = Regen_Energy_ND_Percent*3.6
    
    Input_Energy_ND = InsH2_NP_JCB05_ND(2)
    Input_Energy_ND_Percent = 100*InsH2_NP_JCB05_ND(2)/sum(InsL_JCB05_ND)
    Input_Energy_ND_Pi_Angle = Input_Energy_ND_Percent*3.6
    
    HECM_ND = OutH2_NP_JCB05_ND(3)
    HECM_ND_Percentage = 100*HECM_ND/sum(InsL_JCB05_ND)
    
    Pump_Energy_ND = OutH2_NP_JCB05_ND(2)
    Pump_Energy_ND_Percent = 100*OutH2_NP_JCB05_ND(2)/sum(InsL_JCB05_ND)
    
    Pos_Work_ND = OutH2_NP_JCB05_ND(1)
    Pos_Work_ND_Percent = 100*OutH2_NP_JCB05_ND(1)/sum(InsL_JCB05_ND)
    
    %% For 3 CPR - C0P - Grading 
    Energy_Saving_GR = sum(InsL_JCB05_GR) - sum(InsH3_JCB05_GR)
    Energy_Saving_GR_Percent = 100*Energy_Saving_GR/sum(InsL_JCB05_GR)
    
    Regen_Energy_GR = InsH3_JCB05_GR(1)
    Regen_Energy_GR_Percent = 100*InsH3_JCB05_GR(1)/sum(InsL_JCB05_GR)
    Regen_Energy_GR_Pi_Angle = Regen_Energy_GR_Percent*3.6
    
    Input_Energy_GR = InsH3_JCB05_GR(2)
    Input_Energy_GR_Percent = 100*InsH3_JCB05_GR(2)/sum(InsL_JCB05_GR)
    Input_Energy_GR_Pi_Angle = Input_Energy_GR_Percent*3.6
    
    HECM_GR = OutH3_JCB05_GR(3)
    HECM_GR_Percentage = 100*HECM_GR/sum(InsL_JCB05_GR)
    
    Pump_Energy_GR = OutH3_JCB05_GR(2)
    Pump_Energy_GR_Percent = 100*OutH3_JCB05_GR(2)/sum(InsL_JCB05_GR)
    
    Pos_Work_GR = OutH3_JCB05_GR(1)
    Pos_Work_GR_Percent = 100*OutH3_JCB05_GR(1)/sum(InsL_JCB05_GR)
    
    %% For 2 CPR - C0 - Trenching
    Energy_Saving_TR = sum(InsL_JCB05_TR) - sum(InsH2_NP_JCB05_TR)
    Energy_Saving_TR_Percent = 100*Energy_Saving_TR/sum(InsL_JCB05_TR)
    
    Regen_Energy_TR = InsH2_NP_JCB05_TR(1)
    Regen_Energy_TR_Percent = 100*InsH2_NP_JCB05_TR(1)/sum(InsL_JCB05_TR)
    Regen_Energy_TR_Pi_Angle = Regen_Energy_TR_Percent*3.6
    
    Input_Energy_TR = InsH2_NP_JCB05_TR(2)
    Input_Energy_TR_Percent = 100*InsH2_NP_JCB05_TR(2)/sum(InsL_JCB05_TR)
    Input_Energy_TR_Pi_Angle = Input_Energy_TR_Percent*3.6
    
    HECM_TR = OutH2_NP_JCB05_TR(3)
    HECM_TR_Percentage = 100*HECM_TR/sum(InsL_JCB05_TR)
    
    Pump_Energy_TR = OutH2_NP_JCB05_TR(2)
    Pump_Energy_TR_Percent = 100*OutH2_NP_JCB05_TR(2)/sum(InsL_JCB05_TR)
    
    Pos_Work_TR = OutH2_NP_JCB05_TR(1)
    Pos_Work_TR_Percent = 100*OutH2_NP_JCB05_TR(1)/sum(InsL_JCB05_TR)
    
    %% 3 CPR Unconstraned Sizes No Main Pump - Grading
    Energy_Saving_GR = sum(InsL_JCB05_GR) - sum(InsH3_NP_JCB05_GR)
    Energy_Saving_GR_Percent = 100*Energy_Saving_GR/sum(InsL_JCB05_GR)
    
    Regen_Energy_GR = InsH3_NP_JCB05_GR(1)
    Regen_Energy_GR_Percent = 100*InsH3_NP_JCB05_GR(1)/sum(InsL_JCB05_GR)
    Regen_Energy_GR_Pi_Angle = Regen_Energy_GR_Percent*3.6
    
    Input_Energy_GR = InsH3_NP_JCB05_GR(2)
    Input_Energy_GR_Percent = 100*InsH3_NP_JCB05_GR(2)/sum(InsL_JCB05_GR)
    Input_Energy_GR_Pi_Angle = Input_Energy_GR_Percent*3.6
    
    HECM_GR = OutH3_NP_JCB05_GR(3)
    HECM_GR_Percentage = 100*HECM_GR/sum(InsL_JCB05_GR)
    
    Pump_Energy_GR = OutH3_NP_JCB05_GR(2)
    Pump_Energy_GR_Percent = 100*OutH3_NP_JCB05_GR(2)/sum(InsL_JCB05_GR)
    
    Pos_Work_GR = OutH3_NP_JCB05_GR(1)
    Pos_Work_GR_Percent = 100*OutH3_NP_JCB05_GR(1)/sum(InsL_JCB05_GR)
    