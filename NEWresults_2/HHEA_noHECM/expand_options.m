function [J,lam]=expand_options(options,PR)
% options='C0PM' or 'NPPP' or 'CP00' etc.

%J = ['HECMLosses_c+lam(1)*battery_power_c']; lam=0.15; lamk=1; % WITH battery constraint
J = ['HECMLosses_c']; lam=[]; lamk=0; % NO battery constraint
for k=2:length(PR+1)
    switch options(k)
        case '0'
            lamk=lamk+1;
            J=[J,'+lam(',num2str(lamk),')*QR_c{',num2str(k),'}*PR(',num2str(k),')'];
            lam=[lam,0.1];
        case 'M'
            J = [J,'- QR_c{',num2str(k),'}*PR(',num2str(k),')*R_M_loss(',num2str(k-1),')'];
        case 'P'
            J = [J,'+ QR_c{',num2str(k),'}*PR(',num2str(k),')*R_P_loss(',num2str(k-1),')'];
    end
end
