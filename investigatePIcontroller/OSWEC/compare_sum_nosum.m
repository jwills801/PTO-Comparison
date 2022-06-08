% compare output of reg/irreg test cases both with/without PI using
% a difference in velocity (F_pto = PI(u_actual-u_ideal) or F_pto=PI(u_actual)

x = [150 200];

% load output_reg_nosum/wec_PI_out_A.mat;
% outputAn = output;
% clear output
% 
% load output_reg_sum/wec_PI_out_A.mat;
% outputAs = output;
% clear output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load output_irreg_nosum/wec_PI_out_B.mat;
outputAn = output;
clear output

load output_irreg_sum/wec_PI_out_B.mat;
outputAs = output;
clear output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot(outputAn.controller.time,outputAn.controller.force,'b',...
    outputAs.controller.time,outputAs.controller.force,'r');
legend('no sum','sum');
title('force');
xlim(x);

figure()
plot(outputAn.controller.time,outputAn.controller.velocity,'b',...
    outputAs.controller.time,outputAs.controller.velocity,'r');
legend('no sum','sum');
title('velocity');
xlim(x);

b = 20;
dt = outputAn.controller.time(2)-outputAn.controller.time(1);
ip_an = outputAn.controller.force.*outputAn.controller.velocity;
avep_an = movmean(ip_an,b/dt);
ip_as = outputAs.controller.force.*outputAs.controller.velocity;
avep_as = movmean(ip_as,b/dt);

figure()
plot(outputAn.controller.time,ip_an,'b',...
    outputAn.controller.time,avep_an,'b--',...
    outputAs.controller.time,ip_as,'r',...
    outputAs.controller.time,avep_as,'r--');
legend('no sum instant','no sum average','sum instant','sum average');
title('power');
xlim(x);

figure()
plot(outputAn.controller.time,cumsum(ip_an),'b',...
    outputAs.controller.time,cumsum(ip_as),'r');
legend('Work sum','Work no sum');
title('Work comparison');
% xlim(x);




