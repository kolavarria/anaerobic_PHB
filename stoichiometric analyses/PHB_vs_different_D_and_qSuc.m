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
model = changeRxnBounds(model,'EX_sucrose[e]',-10*(Mw/1000),'l');
model = changeRxnBounds(model,'EX_glc[e]',0*(Mw/1000),'l');
model = changeRxnBounds(model,'EX_o2[e]',0,'l');
model = changeRxnBounds(model,'ACK',0,'b');
%model = changeRxnBounds(model,'PFL',0,'b');
model = changeRxnBounds(model,'EX_etoh[e]',0,'b');
model = changeRxnBounds(model,'EX_lac_D[e]',0,'b');
model = changeRxnBounds(model,'EX_co2[e]',0,'l');
% model = changeRxnBounds(model,'EX_for[e]',0,'b');
model = changeRxnBounds(model,'ATPM',3.2*(Mw/1000),'l');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = changeObjective(model,'ATPM');
%mutante = changeObjective(mutante,'EX_Biomass');

cont=1;
for i=1:24
model = changeRxnBounds(model,'EX_Biomass',i/100,'b');

for j=1:4
    
model = changeRxnBounds(model,'EX_sucrose[e]',-(j+5)*(Mw/1000),'l');    

FBAsolution = optimizeCbModel(model,'max');

if FBAsolution.f
        qSucrose(cont,j)=(-1000/Mw)*FBAsolution.x(findRxnIDs(model,'EX_sucrose[e]'));  
        qATP(cont,j)=FBAsolution.x(findRxnIDs(model,'ATPM'))*1000/Mw;
        qX(cont,j)=FBAsolution.x(findRxnIDs(model,'EX_Biomass'));
        AARflux(cont,j)=(1000/Mw)*FBAsolution.x(findRxnIDs(model,'sintPHB'));
        grams_of_PHB(cont,j)=(86/Mw)*FBAsolution.x(findRxnIDs(model,'sintPHB'));
        PHB_in_total_mass(cont,j)=100*grams_of_PHB(cont,j)/(grams_of_PHB(cont,j)+qX(cont,j));
        PHB_per_lean_biomass(cont,j)=grams_of_PHB(cont,j)/qX(cont,j);
end

end
cont=cont+1;
end

%printFluxVector(model, FBAsolution.x,1,1);

figure(1)
surf(qX,qSucrose,PHB_in_total_mass);
xlabel('D (h^{-1})');
ylabel('q_{sucrose} (mmol/g^{CDW}/h)');
zlabel('PHB content (% of weight)');
colorbar('northoutside');








%figure(2)
% surf(qX,qSucrose,qATP);
% xlabel('D (h^{-1})');
% ylabel('q_{sucrose} (mmol/g^{CDW}/h)');
% zlabel('q_{ATP} (mmol/g^{CDW}/h)');
% colorbar('northoutside');

%figure(3)
% surf(qX,qSucrose,AARflux);
% xlabel('D (h^{-1})');
% ylabel('q_{sucrose} (mmol/g^{CDW}/h)');
% zlabel('flux_{AAR} (mmol/g^{CDW}/h)');
% colorbar('northoutside');
