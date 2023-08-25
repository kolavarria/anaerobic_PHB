clearvars
clc

% Using an external function developed by Stephan Koehler (2022). 
% tsvread: importing tab-separated data (https://www.mathworks.com/matlabcentral/fileexchange/32782-tsvread-importing-tab-separated-data), MATLAB Central File Exchange. 
% Retrieved August 8, 2022.


[data,~, raw] = tsvread('sbtab_ED.tsv');
pathway_name={'ED'};

network=cellfun(@num2str,raw,'un',0); %this step makes easy the reading of the content of the cell

network = table(raw(3:end,1),raw(3:end,2),raw(3:end,3),raw(3:end,8),'VariableNames',{'Parameters' 'Rxns' 'Mets','Values'});

network = table2cell(network); 

[parameterID,parameter_position1,parameterIndex] = unique(network(:,1),'stable');
[rxnID,rxn_position1,rxnIndex] = unique(network(:,2),'stable');

rxnID = rxnID(1:end-1,1);
rxnNames={};
kcatr={};
kcatf={};
KM={};
j1=0;
j2=0;
j3=0;
j4=0;
chain2=[];
record=[];

for i=1:length(parameterIndex)
  
  if parameterIndex(i)==find(contains(parameterID,'equilibrium constant'))
  j1 = j1+1;    
  rxnNames(j1) = network(i,2);   
  end
  
  if parameterIndex(i)==find(contains(parameterID,'product catalytic rate constant'))
  j2 = j2+1;    
  kcatr(j2) = network(i,4);   
  end
  
  if parameterIndex(i)==find(contains(parameterID,'substrate catalytic rate constant'))
  j3 = j3+1;    
  kcatf(j3) = network(i,4);   
  end  
  
  if parameterIndex(i)==find(contains(parameterID,'Michaelis constant'))
      
    if isequal(network(i,2),record)
        chain=[];
        chain=strcat(network(i,3),':',network(i,4));
        chain2 = strcat(chain2,{' '},chain);
        KM{j4}=chain2;
        
    else
        j4=j4+1; record=[]; chain=[]; chain2=[];
        record=network(i,2);
        chain=strcat(network(i,3),':',network(i,4));
        chain2 = strcat(chain2,{' '},chain);
        KM{j4}=chain2;
        
    end
      
  end
  
  
 
end


rxnNames = rxnNames';
kcatr = kcatr';
kcatf = kcatf';
KM = KM';
T1 = table(rxnNames,kcatf,kcatr,KM,'VariableNames',{'Reaction_name' 'kcatf' 'kcatr','KM'});

if isequal(pathway_name,{'EMP'})
D = {'GluKin';'PGI';'PFK';'ALD';'TPI';'GAPDH';'PGK';'PGM';'ENO';'PYK';'AcdDH';'AlcDH';'PFL';'Thio';'AAR'}; % the desired order
end

if isequal(pathway_name,{'ED'})
D = {'GluKin';'G6PDH';'PGase';'GLNdeh';'KDPGald';'GAPDH';'PGK';'PGM';'ENO';'PYK';'AcdDH';'AlcDH';'PFL';'Thio';'AAR'}; % the desired order    
end

if isequal(pathway_name,{'NOGEMP'})
D = {'GluKin';'PGI';'PKT';'TALA';'TKT1';'RibIso';'RibEpi';'TKT2';'TPI';'ALD';'FBPase';'GAPDH';'PGK';'PGM';'ENO';'PYK';'PDH';'PTA';'Thio';'AAR'}; % the desired order
end 

if isequal(pathway_name,{'NOGEMP_suc'})
D = {'SucP';'HPM';'FruKin';'PGI';'PKT';'TALA';'TKT1';'RibIso';'RibEpi';'TKT2';'TPI';'ALD';'FBPase';'GAPDH';'PGK';'PGM';'ENO';'PYK';'PDH';'PTA';'Thio';'AAR'}; % the desired order
end

if isequal(pathway_name,{'NOGf'})
D = {'GluKin';'PGI';'PKT';'TALA';'TKT1';'RibIso';'RibEpi';'TKT2';'TPI';'ALD';'FBPase';'GAPDH';'PGK';'PGM';'ENO';'PYK';'PDH';'PTA';'Thio';'AARf';'PHBs'}; % the desired order
end 

[found, idx] = ismember(T1.Reaction_name, D);
[~, sortorder] = sort(idx);
T1 = T1(sortorder,:);
name=strcat(pathway_name,'.txt');
writetable(T1,string(name),'Delimiter','tab','WriteVariableNames',1);

