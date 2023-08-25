clc
clearvars
load cellular_model

model = addExchangeRxn(model,{'sucrose[e]'},-10,0);
model = addReaction(model,'SucSimporter',{'sucrose[e]','h[e]','sucrose[c]','h[c]'},[-1 -1 1 1],true,0,30);
model = addReaction(model,'SucP',{'sucrose[c]','pi[c]','fru[c]','g1p[c]'},[-1 -1 1 1],true,-30,30);
model = addReaction(model,'FruKin',{'fru[c]','atp[c]','f6p[c]','adp[c]','h[c]'},[-1 -1 1 1 1],true,0,30);
model = addReaction(model,'PenPKT',{'xu5p_D[c]','pi[c]','g3p[c]','actp[c]'},[-1 -1 1 1],true,0,0);
model = addReaction(model,'sintPHB',{'accoa[c]','nadh[c]','h[c]','HB[c]','coa[c]','nad[c]'},[-2 -1 -1 1 2 1],true,0,100);
model = addExchangeRxn(model,{'HB[c]'},0,60);

%VERY IMPORTANT! The constrainst should be introduced in units of mol/Cmol_of_Biomass/h.
% For example, q_glucose of 14 mmol/gCDW/h is represented as -14*(Mw/1000)
model = changeRxnBounds(model,'EX_glc[e]',0*(Mw/1000),'l');
model = changeRxnBounds(model,'ACK',0,'b');
model = changeRxnBounds(model,'EX_etoh[e]',0,'b');
model = changeRxnBounds(model,'EX_lac_D[e]',0,'b');
model = changeRxnBounds(model,'EX_co2[e]',0,'l');
model = changeRxnBounds(model,'EX_Biomass',0,'l');
model = changeRxnBounds(model,'EX_sucrose[e]',-10*(Mw/1000),'l');
model = changeObjective(model,'EX_Biomass');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count_1 = 0;

mATP_1 = 3.2*(Mw/1000);
model = changeRxnBounds(model,'ATPM',mATP_1,'l');

model = changeRxnBounds(model,'EX_o2[e]',0*(Mw/1000),'l');
% FBAsolution = optimizeCbModel(model,'max');
% printFluxVector(model, FBAsolution.x,1,1);

for i=0:5:200
model = changeRxnBounds(model,'EX_o2[e]',(-i/10)*(Mw/1000),'l');

FBAsolution = optimizeCbModel(model,'max');

if FBAsolution.f
        count_1=count_1+1;
        qSucrose_1(count_1)=(-1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_sucrose[e]'));  
        qSucc_1(count_1) = (1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_succ[c]')); 
        growth_rates_1(count_1)=FBAsolution.x(findRxnIDs(model,'EX_Biomass'));
        qH2(count_1) = (1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_H2[c]')); 
        grams_of_PHB_1(count_1)=(86/Mw)*FBAsolution.x(findRxnIDs(model,'sintPHB'));
        percentPHB_1(count_1)=100*grams_of_PHB_1(count_1)/(grams_of_PHB_1(count_1)+growth_rates_1(count_1));
        qO2_1(count_1)=(-1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_o2[e]'));
        mass_yield_1(count_1)=(FBAsolution.x(findRxnIDs(model,'sintPHB'))*86)/(-342.3*FBAsolution.x(findRxnIDs(model,'EX_sucrose[e]')));
end


end

qSucrose_1=qSucrose_1';
mass_yield_1=mass_yield_1';
percentPHB_1=percentPHB_1';
qO2_1=qO2_1';
qSucc_1=qSucc_1';
growth_rates_1=growth_rates_1';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count_2 = 0;
mATP_2 = 16.4*(Mw/1000);
model = changeRxnBounds(model,'ATPM',mATP_2,'l');
for i = 0:5:200
model = changeRxnBounds(model,'EX_o2[e]',(-i/10)*(Mw/1000),'l');

FBAsolution = optimizeCbModel(model,'max');

if FBAsolution.f
        count_2=count_2+1;
        qSucrose_2(count_2)=(-1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_sucrose[e]'));  
        qSucc_2(count_2) = (1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_succ[c]')); 
        growth_rates_2(count_2)=FBAsolution.x(findRxnIDs(model,'EX_Biomass'));
        grams_of_PHB_2(count_2)=(86/Mw)*FBAsolution.x(findRxnIDs(model,'sintPHB'));
        percentPHB_2(count_2)=100*grams_of_PHB_2(count_2)/(grams_of_PHB_2(count_2)+growth_rates_2(count_2));
        qO2_2(count_2)=(-1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_o2[e]'));
        mass_yield_2(count_2)=(FBAsolution.x(findRxnIDs(model,'sintPHB'))*86)/(-342.3*FBAsolution.x(findRxnIDs(model,'EX_sucrose[e]')));
end


end


qSucrose_2=qSucrose_2';
mass_yield_2=mass_yield_2';
percentPHB_2=percentPHB_2';
qO2_2=qO2_2';
qSucc_2=qSucc_2';
growth_rates_2=growth_rates_2';


figure(1)
% graph_title_1=strcat('growth in batch',{' '},'upper row: m_{ATP} =',{' '},num2str(mATP_1*1000/Mw),{' '},'mmol^{ATP}*g_{CDW}*h^{-1}',{' '},'lower row: m_{ATP} =',{' '},num2str(mATP_2*1000/Mw),{' '},'mmol^{ATP}*g_{CDW}*h^{-1}');
% sgt=sgtitle(graph_title_1);
% sgt.FontSize = 8;

subplot(2,2,1);
scatter (qO2_1,growth_rates_1,'k','*');
xlabel('q_{O_2} (mmol*g_{CDW}*h^{-1})');
ylabel('\mu(h^{-1})');
x1 = qO2_1;
y1 = growth_rates_1;
x2 = qO2_1;
y2 = percentPHB_1;
ax1 = gca;
ax1.FontSize = 12;
set(gcf,'color','w');
set(ax1,'XColor','k','YColor','k'); 
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','none','YColor','k');

ax2.FontSize = 12;
hold on 
scatter (qO2_1,percentPHB_1,'k');
xlabel('q_{O_2} (mmol*g_{CDW}*h^{-1})');
ylabel('PHB content (% of total biomass)');

subplot(2,2,2);

scatter (qO2_1,qSucc_1,'k','x');
ylabel('q_{succinate} (mmol*g_{CDW}*h^{-1})')
xlabel('q_{O_2} (mmol*g_{CDW}*h^{-1})');
x1 = qO2_1;
y1 = qSucc_1;
x2 = qO2_1;
y2 = mass_yield_1;
ax1 = gca;
ax1.FontSize = 12;
set(gcf,'color','w');
set(ax1,'XColor','k','YColor','k'); 
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','none','YColor','k');

ax2.FontSize = 12;
hold on 
scatter (qO2_1,mass_yield_1,'k','+');
ylabel('yield (g^{PHB}/g^{sucrose})');


subplot(2,2,3);
scatter (qO2_2,growth_rates_2,'k','*');
xlabel('q_{O_2} (mmol*g_{CDW}*h^{-1})');
ylabel('\mu(h^{-1})');
x1 = qO2_2;
y1 = growth_rates_2;
x2 = qO2_2;
y2 = percentPHB_2;
ax1 = gca;
ax1.FontSize = 12;
set(gcf,'color','w');
set(ax1,'XColor','k','YColor','k'); 
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','none','YColor','k');

ax2.FontSize = 12;
hold on 
scatter (qO2_2,percentPHB_2,'k');
xlabel('q_{O_2} (mmol*g_{CDW}*h^{-1})');
ylabel('PHB content (% of total biomass)');

subplot(2,2,4);

scatter (qO2_2,qSucc_2,'k','x');
ylabel('q_{succinate} (mmol*g_{CDW}*h^{-1})')
xlabel('q_{O_2} (mmol*g_{CDW}*h^{-1})');
x1 = qO2_2;
y1 = qSucc_2;
x2 = qO2_2;
y2 = mass_yield_2;
ax1 = gca;
ax1.FontSize = 12;
set(gcf,'color','w');
set(ax1,'XColor','k','YColor','k'); 
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','none','YColor','k');

ax2.FontSize = 12;
hold on 
scatter (qO2_2,mass_yield_2,'k','+');
ylabel('yield (g^{PHB}/g^{sucrose})');
