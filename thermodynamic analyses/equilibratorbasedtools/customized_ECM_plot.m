clearvars
clc

rxns = readtable('molarNOGEMP.csv');

costs_table = readtable('costsNOGEMP.csv');

costs=table2array(costs_table);

rxns_names=table2cell(rxns(:,1));

%load data_enzyme_cost

%2D array with basic capacity concentrations (blue bars) on 1st column, fold increase in demand due to thermodynamics (orange bar) in the 2nd column and fold increase in demand due to saturation (brown bar) in 3rd. 
%The 4th column is related to regulation fold change, which was not analyzed in our study.

reactions_CGP = categorical(string(rxns_names));
reactions_CGP = reordercats(reactions_CGP,rxns_names);

[rows,columns]=size(costs);

all_costs=zeros(rows,3);


for i=1:rows
    accumulated=1;
    for j=1:columns-1
        cost_rxn(i)=accumulated*costs(i,j);
        accumulated = cost_rxn(i);
        all_costs(i,j) = cost_rxn(i);
    end
    
end

cost_rxn=cost_rxn';

bh=bar(reactions_CGP,all_costs,'stacked');
set(bh, 'FaceColor', 'Flat')
bh(1).CData = [0 0.4470 0.7410];  % Change color to first level
bh(2).CData = [0.9290 0.6940 0.1250];  % Change color to second level
bh(3).CData = [0.8500 0.3250 0.0980];  % Change color to third level

set(gca, 'YScale', 'log')
ylim([10^-9 10^-2]);
ylabel('required enzyme (M)')
legend('turnover','reversibility','saturation','Location','northwest');
set(gcf,'color','w');
