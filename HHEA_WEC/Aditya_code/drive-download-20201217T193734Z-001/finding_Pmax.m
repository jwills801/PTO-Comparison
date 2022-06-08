load JCB_LongerDC.mat
load 28point1ccnolock_A6VM_4Q_20191210.mat


% cycle = NinetyDeg;
% cycle = Grading;
cycle = Trenching;

D = Disp_28cc;
D_HECM4 = 18.1*(1/100)^3 ;

ScaleHECM4 = D_HECM4/(D*2*pi);

maxP = 0;

lrp = 0.08; %percentage of lower rail pressure wrt to Pmax (Pmax to be found)

spacing = 1;

t_all = 1:spacing:length(cycle.time);
this_chunk = t_all;
F1 = -cycle.boom_f(this_chunk);
ACap1 = boom_Ac;
ARod1 = boom_Ar;

F2 = -cycle.arm_f(this_chunk);
ACap2 = arm_Ac;
ARod2 = arm_Ar;

F3 = -cycle.bkt_f(this_chunk);
ACap3 = bkt_Ac;
ARod3 = bkt_Ar;

T4 = -cycle.sw_t(this_chunk)/Swing_ratio;

P = max(F1)/(ARod1 - lrp*ACap1);
% P = max(F1)/(ARod1);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

P = min(F1)/(lrp*ARod1 - ACap1);
% P = min(F1)/(-ACap1);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

P = max(F2)/(ARod2 - lrp*ACap2);
% P = max(F2)/(ARod2);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

P = min(F2)/(lrp*ARod2 - ACap2);
% P = min(F2)/(ACap2);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

P = max(F3)/(ARod3 - lrp*ACap3);
% P = max(F3)/(ARod3);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

P = min(F3)/(lrp*ARod3 - ACap3);
% P = min(F3)/(ACap3);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

P = max(T4)/(D*ScaleHECM4);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

P = min(T4)/(D*ScaleHECM4);
disp(P)
if abs(P) > maxP
    maxP = abs(P);
end

maxP
    
   
    
    
    