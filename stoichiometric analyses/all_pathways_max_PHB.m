% Author: Karel Olavarria Gamez

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
%Loading the glycolytic model
load glycolytic_model

model = addExchangeRxn(model,{'Suc'},-30,0);
model = addReaction(model,'SucP',{'Suc','Pi','G1P','Fruc'},[-1 -1 1 1],true,0,30);
model = addReaction(model,'HPM',{'G1P','G6P'},[-1 1],true,-30,30);
model = addReaction(model,'FruKin',{'Fruc','ATP','F6P','ADP'},[-1 -1 1 1],true,0,30);

%Initial state
model = changeRxnBounds(model,'EX_Gluc',-1,'l');
model = changeRxnBounds(model,'EX_Suc',0,'b');
model = changeRxnBounds(model,'PFK',0,'l');
model = changeRxnBounds(model,'G6PDH',0,'l');
model = changeRxnBounds(model,'PKT',0,'b');
model = changeRxnBounds(model,'PDH',0,'l');
model = changeRxnBounds(model,'PYK',0,'l');
model = changeRxnBounds(model,'PGase',0,'l');
model = changeRxnBounds(model,'GND',0,'l');
model = changeRxnBounds(model,'FBPase',0,'l');
model = changeRxnBounds(model,'ATPM',0,'l');
model = changeObjective(model,'EX_PHB');

for j=1:5
    
    %EMP
    if j==1
        pathway={'EMP'};
        model = changeRxnBounds(model,'G6PDH',0,'b');
        model = changeRxnBounds(model,'EX_Suc',0,'b');
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
    
% Optimization

 FBAsolution = optimizeCbModel(model,'max');
 
%Flux visualization  
    
  if j==1  
      
        load map_EMP
        rxnID = findRxnIDs(model,reactions);

        for i=1:length(rxnID)    
            fluxes(i)=FBAsolution.x(rxnID(i));
            labels{i} = num2str(sprintf('%g',round(fluxes(i)*100)/100));
            text(axis_x(i),axis_y(i),labels(i),'FontSize',8,'Color','red','FontWeight','bold');
        end

  end       
        
  if j==2  
      
        load map_ED
        rxnID = findRxnIDs(model,reactions);

        for i=1:length(rxnID)    
            fluxes(i)=FBAsolution.x(rxnID(i));
            labels{i} = num2str(sprintf('%g',round(fluxes(i)*100)/100));
            text(axis_x(i),axis_y(i),labels(i),'FontSize',8,'Color','red','FontWeight','bold');
        end

  end            
  
  if j==3  
      
        load map_NOGEMP
        rxnID = findRxnIDs(model,reactions);

        for i=1:length(rxnID)    
            fluxes(i)=FBAsolution.x(rxnID(i));
            labels{i} = num2str(sprintf('%g',round(fluxes(i)*100)/100));
            text(axis_x(i),axis_y(i),labels(i),'FontSize',8,'Color','red','FontWeight','bold');
        end

  end  
  
   if j==4  
      
        load map_NOGEMPsuc
        rxnID = findRxnIDs(model,reactions);

        for i=1:length(rxnID)    
            fluxes(i)=FBAsolution.x(rxnID(i));
            labels{i} = num2str(sprintf('%g',round(fluxes(i)*100)/100));
            text(axis_x(i),axis_y(i),labels(i),'FontSize',8,'Color','black','FontWeight','bold');
        end

  end
  
        hold on
        
clear reactions rxnIDs fluxes labels
end
