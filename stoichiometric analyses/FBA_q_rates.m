% Author: Karel Olavarria Gamez

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
%Loading the glycolytic model
load glycolytic_model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = addReaction(model,'SucP',{'Suc','Pi','G1P','Fruc'},[-1 -1 1 1],true,0,30);
model = addReaction(model,'HPM',{'G1P','G6P'},[-1 1],true,-30,30);
model = addReaction(model,'FruKin',{'Fruc','ATP','F6P','ADP'},[-1 -1 1 1],true,0,30);


%Initial state
model = changeRxnBounds(model,'PFK',0,'l');
model = changeRxnBounds(model,'G6PDH',0,'l');
model = changeRxnBounds(model,'PKT',0,'b');
model = changeRxnBounds(model,'PDH',0,'l');
model = changeRxnBounds(model,'PYK',0,'l');
model = changeRxnBounds(model,'PGase',0,'l');
model = changeRxnBounds(model,'GND',0,'l');
model = changeRxnBounds(model,'FBPase',0,'l');
model = changeRxnBounds(model,'ATPM',0,'l');

%selecting Objective
model = changeObjective(model,'PHBs');

%Empty Table to be filled in the loop
global_table = table;

for j=1:4

    %EMP
    if j==1
        pathway={'EMP'};
        model = changeRxnBounds(model,'EX_Suc',0,'b');
        model = changeRxnBounds(model,'EX_Gluc',-1,'l');
        model = changeRxnBounds(model,'G6PDH',0,'b');
    end

    %EDP
    if j==2
         pathway={'ED'};
         model = changeRxnBounds(model,'EX_Suc',0,'b');
         model = changeRxnBounds(model,'EX_Gluc',-1,'l');
         model = changeRxnBounds(model,'PFK',0,'b');
         model = changeRxnBounds(model,'G6PDH',0,'l'); 
         model = changeRxnBounds(model,'G6PDH',100,'u');  
    end

    %NOG
    if j==3
         pathway={'NOGEMP'};
         model = changeRxnBounds(model,'EX_Suc',0,'b');
         model = changeRxnBounds(model,'EX_Gluc',-1,'l');
         model = changeRxnBounds(model,'PKT',0,'l');
         model = changeRxnBounds(model,'PKT',100,'u');
         model = changeRxnBounds(model,'GluKin',0,'l');
         model = changeRxnBounds(model,'GluKin',100,'u');
         model = changeRxnBounds(model,'G6PDH',0,'b');
         model = changeRxnBounds(model,'PFK',0,'b');
         model = changeRxnBounds(model,'AlcDH',0,'b');
         model = changeRxnBounds(model,'PFL',0,'b');
         model = removeRxns(model,{'ALD','TPI'});
         model = addReaction(model,'ALD','DHAP + G3P <=> FBP');
         model = addReaction(model,'TPI','G3P <=> DHAP');
    end
    
    %NOGsuc
    if j==4
         pathway={'NOGEMPsuc'};
         model = changeRxnBounds(model,'EX_Suc',-1,'l');
         model = changeRxnBounds(model,'EX_Gluc',0,'b');
         model = changeRxnBounds(model,'PKT',0,'l');
         model = changeRxnBounds(model,'PKT',100,'u');
         model = changeRxnBounds(model,'GluKin',0,'l');
         model = changeRxnBounds(model,'GluKin',100,'u');
         model = changeRxnBounds(model,'G6PDH',0,'b');
         model = changeRxnBounds(model,'PFK',0,'b');
         model = changeRxnBounds(model,'AlcDH',0,'b');
         model = changeRxnBounds(model,'PFL',0,'b');
    end
    
            % Optimizing 
            FBAsolution = optimizeCbModel(model,'max');
            disp ('-----------------------------')
            title = strcat('Global Balance for the mode_',string(pathway));
            disp(title)
            disp ('-----------------------------')
            printFluxVector(model, FBAsolution.x,1,1);

    %Generating the file containing the q-rates
    qGluc = round(FBAsolution.x(findRxnIDs(model,'EX_Gluc')) * 10)/10;
    qEt = round(FBAsolution.x(findRxnIDs(model,'EX_ETOH')) * 10)/10;
    qFor = round(FBAsolution.x(findRxnIDs(model,'EX_FOR')) * 10)/10;
    qATP = round( abs(FBAsolution.x(findRxnIDs(model,'ATPM')) ) * 10)/10;
    qPHB = round(FBAsolution.x(findRxnIDs(model,'EX_PHB')) * 10)/10;
    qPHB_theo_max = round(( abs(qGluc*24)/36) * 10)/10;
    percent_of_Y_theor = round(((qPHB/qPHB_theo_max)*100) * 10)/10;
    T2 = table(qGluc',qEt',qFor',qATP',qPHB',percent_of_Y_theor','VariableNames',{'qGluc' 'qEt' 'qFor' 'qATP' 'qPHB' 'Percent of Y_max'});
    global_table = vertcat(global_table,T2);
       
clear T2
end

writetable(global_table,'q_rates.txt','Delimiter','tab');


disp('-------------------');
disp('Whole set of reactions');
disp('-------------------');

printRxnFormula(model);
reaction_formula=ans;
reaction_names=model.rxns;
T3=table(reaction_names,reaction_formula);
writetable(T3,'all_reactions.txt','Delimiter','tab');

