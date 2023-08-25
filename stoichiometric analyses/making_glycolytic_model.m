% Author: Karel Olavarria Gamez

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
model = createModel();

% Taking up the carbon source

model = addReaction(model,'GluKin',{'Gluc','ATP','G6P','ADP'},[-1 -1 1 1],true,0,30);
model = changeRxnBounds(model,'GluKin',0,'l');

% Embden-Meyerhof pathway

model = addReaction(model,'PGI','G6P <=> F6P');
model = addReaction(model,'PFK','ATP + F6P <=> ADP + FBP');
model = changeRxnBounds(model,'PFK',0,'l');
model = addReaction(model,'ALD','FBP <=> DHAP + G3P');
model = addReaction(model,'TPI','DHAP <=> G3P');
model = addReaction(model,'GAPDH','G3P + NAD + Pi <=> BPG + NADH');
model = addReaction(model,'PGK','BPG + ADP <=> P3G + ATP');
model = addReaction(model,'PGM','P3G <=> P2G');
model = addReaction(model,'ENO','P2G <=> H2O + PEP');
model = addReaction(model,'PYK','ADP + PEP <=> ATP + PYR');
model = changeRxnBounds(model,'PYK',0,'l');
model = addReaction(model,'PDH','CoA + NAD + PYR <=> AcCoA + CO2 + NADH');
model = changeRxnBounds(model,'PDH',0,'l');

%Phosphoketolase and Phosphotransacetylase
model = addReaction(model,'PKT','X5P + Pi <=> G3P + AcP + H2O');
model = changeRxnBounds(model,'PKT',0,'b');
model = addReaction(model,'PTA','AcP + CoA <=> AcCoA + Pi');

% Pentose-Phosphate pathway

model = addReaction(model,'G6PDH','G6P + NAD <=> P6Glac + NADH');
model = changeRxnBounds(model,'G6PDH',0,'l');
model = addReaction(model,'PGase','P6Glac + H2O <=> P6GLN');
model = changeRxnBounds(model,'PGase',0,'l');
model = addReaction(model,'GND','P6GLN + NAD <=> CO2 + NADH + Ru5P');
model = changeRxnBounds(model,'GND',0,'l');
model = addReaction(model,'RibIso','Ri5P <=> Ru5P');
model = addReaction(model,'RibEpi','Ru5P <=> X5P');
model = addReaction(model,'TALA','E4P + F6P <=> G3P + S7P');
model = addReaction(model,'TKT1','G3P + S7P <=> Ri5P + X5P');
model = addReaction(model,'TKT2','F6P + G3P <=> E4P + X5P');

%Gluconeogenic reactions
model = addReaction(model,'FBPase','FBP + H2O <=> F6P + Pi');
model = changeRxnBounds(model,'FBPase',0,'l');

% Entner-Doudoroff pathway
model = addReaction(model,'GLNdeh',{'P6GLN','H2O','KDPG'},[-1 1 1],true,0,100);
model = addReaction(model,'KDPGald',{'KDPG','PYR','G3P'},[-1 1 1],true,0,100);

% Fermentation pathways
model = addReaction(model,'AlcDH','AcCoA + 2 NADH <=> CoA + ETOH + 2 NAD');
model = addReaction(model,'PFL','PYR + CoA <=> AcCoA + FOR');

% Reaction enabling NADH-dependent PHB accumulation
model = addReaction(model,'Thio','AcCoA + AcCoA <=> AcAcCoA + CoA');
model = addReaction(model,'AAR','AcAcCoA + NADH <=> HBCoA + NAD');
model = addReaction(model,'PHBs','HBCoA + HBCoA + H2O <=> PHB + CoA + CoA');

% ATP hydrolysis
model = addReaction(model,'ATPM','ATP + H2O -> ADP + Pi');

% Exchange Rections
model = addExchangeRxn(model,{'Gluc'},-100,0);
model = addExchangeRxn(model,{'Ace'},0,100);
model = addExchangeRxn(model,{'CO2'},0,100);
model = addExchangeRxn(model,{'FOR'},0,100);
model = addExchangeRxn(model,{'ETOH'},0,100);
model = addExchangeRxn(model,{'LAC'},0,100);
model = addExchangeRxn(model,{'H2O'},-100,100);
model = addExchangeRxn(model,{'Pi'},-100,100);
model = addExchangeRxn(model,{'PHB'},0,100);

maintenance_ATP = 0;
model = changeRxnBounds(model,'ATPM',maintenance_ATP,'l');

%selecting Objective
model = changeObjective(model,'ATPM');

clear maintenance_ATP
save glycolytic_model
model_sbml = writeCbModel(model,'format','sbml');