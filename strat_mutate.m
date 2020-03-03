function [ mut_strat_test ] = strat_mutate(which_sensor,which_gene_mut,mother_strat )
% 
%[fnull | fN  UN | fTB  UTB | fTA  UTA | fQS QS]
mut_strat_test = mother_strat;
genes_possible_mut = [1,2*find(which_sensor),2*find(which_sensor)+1];
gene_to_mut = genes_possible_mut(which_gene_mut);
% perform mutation mutate
mut_strat_test(gene_to_mut) = mother_strat(gene_to_mut) + normrnd(0,0.001);
end
 
