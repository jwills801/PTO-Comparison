% Rail optimization brute force for combined torques for both actuators
% Torque is minimized. 4 rails assumed (2 middle ones to be optimized)

nrails = 4;

% regu;ar waves
%Force limits normalized to Pmax*A
F_limits = [-0.9236, .7036]; %defined at 35MPa
%F_limits = [-1 2; 0 2]; % no limit

Aratio = 1; % area ratio

p = 0:0.01:1;
dF = ones(1,length(p));
N = length(p);
least = 100;

switch nrails
    case 3
        for k1=2:N-1
            ind = k1;
            rails = [0,p(ind),1]';
            maxDF1 = gap_to_go_limit(Aratio(1), rails, F_limits(1,:),1);
            if maxDF1<least
               least = maxDF1; 
               least_ind = ind; 
            end
        end
    case 4
        for k1=2:N-1
            for k2=k1+1:N-1
                ind = [k1,k2];
                rails = [0,p(ind),1]';
                maxDF1 = gap_to_go_limit(Aratio(1), rails, F_limits(1,:),1);
                if maxDF1<least
                   least = maxDF1; 
                   least_ind = ind;
                end
           end
        end
    otherwise
        ME = MExecption('Invalid nrails - must be 3 or 4');
        throw(ME);
end

rails = [0,p(least_ind),1]';
disp('Optimal rails: '); disp(rails')
least = gap_to_go_limit(Aratio(1), rails, F_limits(1,:),1);
rails = linspace(0,1,nrails)';
uniform = gap_to_go_limit(Aratio(1), rails, F_limits(1,:),1);
disp('Optimal and uniform')
disp([least,uniform])




function worst = gap_to_go_limit(Aratio, rails,limits,cavit)
% Aratio = area ratio
% normalized rail pressures as column: [0, p(2), p(3), 1]'
% limits = [min F, max F]
    nrails = length(rails);
    DF1s = Aratio*rails*ones(1,nrails)-ones(nrails,1)*rails';
    DFs = sort(unique(DF1s(:)));
    DF1s = sort(DF1s(:)); 
    a = max(1,find(DFs>limits(1),1,'first')-1);
    b = min(length(DFs),find(DFs<limits(2),1,'last')+1);
    DFs = DFs(a:b);
    DF1s = DF1s(DF1s>=DFs(1) & DF1s<=DFs(end));
    worst = max(diff(DFs))/2;
    if cavit
    cav = 0;
    for j=1:nrails
        if sum(DF1s(:)==Aratio*rails(j))==1 %cavitating rail
            tmp=DFs-Aratio*rails(j); 
            %Gap above is the smallest positive tmp
            %Gap between below force and above force is the |small negative tmp| + Gap above 
            %CODE ALERT: Only valid if rail force below the cavitating
            %rail force is not itself a cavitating rail
            cav_ = [(min(abs(tmp(tmp<0)))+min(tmp(tmp>0)))/2,min(tmp(tmp>0))];
            cav = max([cav,min(cav_)]);
        end
    end
    worst= max(worst, cav);
    end
end
