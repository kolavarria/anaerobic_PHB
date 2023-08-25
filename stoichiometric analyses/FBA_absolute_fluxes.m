% Author: Karel Olavarria Gamez

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
%Loading the glycolytic model
load glycolytic_model

% Taking up the carbon source

model = addReaction(model,'SucP',{'Suc','Pi','G1P','Fruc'},[-1 -1 1 1],true,0,30);
model = addReaction(model,'HPM',{'G1P','G6P'},[-1 1],true,-30,30);
model = addReaction(model,'FruKin',{'Fruc','ATP','F6P','ADP'},[-1 -1 1 1],true,0,30);


% Exchange Rections
model = addExchangeRxn(model,{'Suc'},-10,0);

%Initial state
model = changeRxnBounds(model,'PFK',0,'l');
model = changeRxnBounds(model,'G6PDH',0,'l');
model = changeRxnBounds(model,'PKT',0,'b');
model = changeRxnBounds(model,'PDH',0,'l');
model = changeRxnBounds(model,'PYK',0,'l');
model = changeRxnBounds(model,'PGase',0,'l');
model = changeRxnBounds(model,'GND',0,'l');
model = changeRxnBounds(model,'FBPase',0,'l');
maintenance_ATP = 3.2;
model = changeRxnBounds(model,'ATPM',maintenance_ATP,'l');

%Empty Table to be filled in the loop
global_table = table;

for j=1:4

    %EMP
    if j==1
        pathway={'EMP'};
        model = changeRxnBounds(model,'G6PDH',0,'b');
        model = changeRxnBounds(model,'EX_Suc',0,'b'); 
        model = changeObjective(model,'EX_Gluc');
    end

    %EDP
    if j==2
         pathway={'ED'};
         model = changeRxnBounds(model,'EX_Suc',0,'b'); 
         model = changeRxnBounds(model,'PFK',0,'b');
         model = changeRxnBounds(model,'G6PDH',0,'l'); 
         model = changeRxnBounds(model,'G6PDH',100,'u');  
    end

    %NOG
    if j==3
         pathway={'NOGEMP'};
         model = changeRxnBounds(model,'EX_Suc',0,'b'); 
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
         model = changeRxnBounds(model,'EX_Suc',-10,'l');
         model = changeRxnBounds(model,'EX_Gluc',0,'b');
         model = changeObjective(model,'EX_Suc');
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

   % Writing the output files
    counter = 0; ActiveRxns = {};

            for i=1:length(model.rxns)

                if FBAsolution.x(i)<-0.01 || FBAsolution.x(i)>0.01
                    counter = counter+1;
                    ActiveRxns(counter) = model.rxns(i);
                
                        if j==4
                            Relative_Fluxes(counter) = -1*FBAsolution.x(i)/FBAsolution.x(findRxnIDs(model,'EX_Suc'));
                        else
                            Relative_Fluxes(counter) = -1*FBAsolution.x(i)/FBAsolution.x(findRxnIDs(model,'EX_Gluc'));    
                        end
                    
                    Absolute_Fluxes(counter) = FBAsolution.x(i)*(1/1000)*(1/0.0019)*(1/3600);
                end

            end

    ActiveRxnsFormula = printRxnFormula(model,ActiveRxns,false); %clc;

    %Generating the file containing the reactions carrying fluxes
    T1 = table(ActiveRxns',ActiveRxnsFormula',Relative_Fluxes',Absolute_Fluxes','VariableNames',{'Reaction Name' 'Reaction Formula' 'Relative Flux' 'Absolute Flux (M/s)'});
    writetable(T1,strcat('Active Fluxes_',string(pathway),'.txt'),'Delimiter','tab');
      
clear ActiveRxnsFormula Relative_Fluxes Absolute_Fluxes ActiveRxns ActiveRxnsFormula T1 T2 reactions rxnIDs fluxes labels
end