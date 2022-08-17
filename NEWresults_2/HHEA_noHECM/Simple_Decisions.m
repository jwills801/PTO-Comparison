for i = 1:length(t)-1
    Decision_vector(i) = find(abs ( Frange(:) - (-F1(i+1)) ) < 1e-3 ,1);
end

figure, plot(t(1:end-1),Decision_vector)
starting_ind = Decision_vector(1);

d_ind = NaN(1,length(t));
d_ind(1) = sub2ind([length(t),length(PRA1)],1,Decision_vector(1));
for t_ind = 2:length(t)
    d_ind(t_ind) = sub2ind([length(t),length(PRA1)],t_ind,Decision_vector(t_ind-1));
end

Decision_vector_c = Decision_vector;
d_ind_c = d_ind;