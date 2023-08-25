% Author: Karel Olavarria Gamez
clc
clearvars
%Loading the glycolytic model
load glycolytic_model
maintenance_ATP = 0;
model = changeRxnBounds(model,'ATPM',maintenance_ATP,'l');

%selecting Objective
model = changeObjective(model,'ATPM');

%Empty Table to be filled in the loop
global_table = table;


load map_EMPED
s = warning('off', 'Images:initSize:adjustingMag');
imshow(img, 'InitialMagnification',100, 'Border','tight')
warning(s);
figureHandle = gcf;

% To print fluxes on a metabolic map
rxnID = findRxnIDs(model,reactions);
colors=[[1 0 0]; [0 0 1]];

for j=1:2

    %EMP
    if j==1
    model = changeRxnBounds(model,'G6PDH',0,'b');     
    model = changeRxnBounds(model,'EX_Gluc',-1,'l');
    end

    %EDP
    if j==2
    model = changeRxnBounds(model,'PFK',0,'b');
    model = changeRxnBounds(model,'G6PDH',0,'l'); 
    model = changeRxnBounds(model,'G6PDH',100,'u');  
    model = changeRxnBounds(model,'EX_Gluc',-1,'l');
    end
   
            % Optimizing ATP production 
            FBAsolution = optimizeCbModel(model,'max');
            disp ('-----------------------------')
            title = strcat('Global Balance for the mode_',num2str(j));
            disp(title)
            disp ('-----------------------------')
            printFluxVector(model, FBAsolution.x,1,1);
            
   % Writing the output files
    counter = 0; ActiveRxns = {};

            for i=1:length(model.rxns)

                if FBAsolution.x(i)<-0.01 || FBAsolution.x(i)>0.01
                counter = counter+1;
                ActiveRxns(counter) = model.rxns(i);
                Relative_Fluxes(counter) = -1*FBAsolution.x(i)/FBAsolution.x(findRxnIDs(model,'EX_Gluc'));
                Absolute_Fluxes(counter) = FBAsolution.x(i)*(1/1000)*(1/0.0019)*(1/3600);
                end

            end

    ActiveRxnsFormula = printRxnFormula(model,ActiveRxns,false); %clc;

    %Generating the file containing the reactions carrying fluxes
    T1 = table(ActiveRxnsFormula',Relative_Fluxes',Absolute_Fluxes',ActiveRxns','VariableNames',{'Reaction Formula' 'Relative Flux' 'Absolute Flux (M/s)' 'Reaction Name'});
    writetable(T1,strcat('Active Fluxes_',num2str(j),'.txt'),'Delimiter','tab');
    %Generating the file containing the q-rates
    qGluc = round(FBAsolution.x(findRxnIDs(model,'EX_Gluc')) * 10)/10;
    qEt = round(FBAsolution.x(findRxnIDs(model,'EX_ETOH')) * 10)/10;
    qFor = round(FBAsolution.x(findRxnIDs(model,'EX_FOR')) * 10)/10;
    qAce = round(FBAsolution.x(findRxnIDs(model,'EX_Ace')) * 10)/10;
    qATP = round( abs(FBAsolution.x(findRxnIDs(model,'ATPM')) ) * 10)/10;
    qPHB = round(FBAsolution.x(findRxnIDs(model,'EX_PHB')) * 10)/10;
    qPHB_theo_max = round(( abs(qGluc*24)/36) * 10)/10;
    percent_of_Y_theor = round(((qPHB/qPHB_theo_max)*100) * 10)/10;
    
    T2 = table(qGluc',qEt',qFor',qATP',qPHB',percent_of_Y_theor','VariableNames',{'qGluc' 'qEt' 'qFor' 'qATP' 'qPHB' 'Percent of Y_max'});
    global_table = vertcat(global_table,T2);
    
        
        for i=1:length(rxnID)    
            fluxes(i)=FBAsolution.x(rxnID(i));
        end

        for i = 1:length(reactions)  
            
            if fluxes(i)==0
                labels{i} = num2str(0);
            else
                labels{i} = num2str(sprintf('%g',round(fluxes(i)*100)/100));
            end
            
        end

        for i=1:length(reactions)
        text(axis_x(i),axis_y(i)+20*(j-1),labels(i),'FontSize',8,'Color',colors(j,:));
        end

        hold on

clear ActiveRxnsFormula Relative_Fluxes Absolute_Fluxes ActiveRxns ActiveRxnsFormula T1 T2
end

writetable(global_table,'Balances.txt','Delimiter','tab');

