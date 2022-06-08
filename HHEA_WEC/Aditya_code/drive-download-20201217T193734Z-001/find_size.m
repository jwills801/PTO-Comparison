% Size components
load JCB_LongerDC.mat

names = {'NinetyDeg','Grading','Trenching'};
% maxRPM = 3000; %RPM
for k=1:3
    eval(['cycle_FS=',names{k},';'])
    D_boom(k)=max(abs(cycle_FS.boom_v))*boom_Ar/(maxRPM/60);
    D_arm(k)=max(abs(cycle_FS.arm_v))*arm_Ar/(maxRPM/60);
    D_bkt(k)=max(abs(cycle_FS.bkt_v))*bkt_Ar/(maxRPM/60);
end
% Displacements in [m^3/rev]
D_HECM1=max(D_boom);
D_HECM2=max(D_arm);
D_HECM3=max(D_bkt);

% TO FIND D_HECM4, 1) find Pmax for each cycle; 2) Find minimum
% displacement for each cycle; 3) Find max of the displacement required for
% all cycles




