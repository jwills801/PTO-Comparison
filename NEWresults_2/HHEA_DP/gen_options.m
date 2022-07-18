function all_options = gen_options(nrails, choices)
% Find all options for given number of CPRs 
% nrails  = number of CPRs. Valid options are [3, 4, 5].
% Choices = "['0','M','P']" normally but can be reduced to not generate unlikely cases. 
% Options with no 'P' are eliminated.

%choices=['0','M','P'];

k=0;
switch nrails
    case 3 %3 rail options
    for k1=1:length(choices)
        for k2=1:length(choices)
            k=k+1;
            if (choices(k1)=='0' | choices(k2) == '0'),
                all_options{k}=['C',choices(k1),choices(k2)];
            else 
                all_options{k}=['N',choices(k1),choices(k2)];
            end
            if ~ismember('P',all_options{k}), k=k-1; end % eliminate the choice.
        end
    end

    case 4  %4 rail options
    for k1=1:length(choices)
        for k2=1:length(choices)
            for k3=1:length(choices)
                k=k+1;
                if (choices(k1)=='0' | choices(k2) == '0' | choices(k3) == '0'),
                    all_options{k}=['C',choices(k1),choices(k2),choices(k3)];
                else 
                    all_options{k}=['N',choices(k1),choices(k2),choices(k3)];
                end
                if ~ismember('P',all_options{k}), k=k-1; end % eliminate the choice.
            end
        end
    end
    case 5
    for k1=1:length(choices)
        for k2=1:length(choices)
            for k3=1:length(choices)
                for k4=1:length(choices)
                    k=k+1;
                    if (choices(k1)=='0' | choices(k2) == '0' | choices(k3) == '0'),
                        all_options{k}=['C',choices(k1),choices(k2),choices(k3),choices(k4)];
                    else 
                        all_options{k}=['N',choices(k1),choices(k2),choices(k3),choices(k4)];
                    end
                    if ~ismember('P',all_options{k}), k=k-1; end % eliminate the choice.
                end
            end
        end
    end
end
end

