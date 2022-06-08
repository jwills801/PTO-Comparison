function [MP,J,lam]=MYexpand_options(options,PR)
% options='C0PM' or 'NPPP' or 'CP00' etc.

J = ['min(MPLosses + HECMLosses_noisy']; % + lam(1)*battery_power'];
MP = ['MPLosses = zeros(size(HECMLosses))'];
lam=[];

lamk=0;
for k=2:length(PR)
    switch options(k)
        case '0'
            lamk=lamk+1;
            J=[J,' + lam(',num2str(lamk),')*QR{',num2str(k),'}']; % *PR(',num2str(k),' That is the flow constraint
            lam=[lam,-.3];
        case 'M'
            MP = [MP,'- QR{',num2str(k),'}*PR(',num2str(k),')*R_M_loss(',num2str(k-1),')'];
        case 'P'
            MP = [MP,'+ QR{',num2str(k),'}*PR(',num2str(k),')*R_P_loss(',num2str(k-1),')'];
    end
end
J = [J, ' , [] , 2)'];
MP = [MP, ';'];

% Protype
% case 'C0PM'
%     MPLosses = QR{3}*PR(3)*R_P_loss(2) - QR{4}*PR(4)*R_M_loss(3);
%     f = @(lam) -sum(min(MPLosses + HECMLosses+lam(1)*battery_power+lam(2)*QR{2}*PR(2),[],2));
%     lam = [0.15,0.1];
