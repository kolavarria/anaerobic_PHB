%Author: Karel Olavarria Gamez
clc
clearvars
model = createModel();

model = addExchangeRxn(model,{'so4[e]'},-200,200);
model = addExchangeRxn(model,{'ac[e]'},0,200);
model = addExchangeRxn(model,{'co2[e]'},0,200);
model = addExchangeRxn(model,{'etoh[e]'},0,200);
model = addExchangeRxn(model,{'for[e]'},0,200);
model = addExchangeRxn(model,{'glc[e]'},-200,0);
model = addExchangeRxn(model,{'h[e]'},-200,200);
model = addExchangeRxn(model,{'h2o[e]'},-200,200);
model = addExchangeRxn(model,{'lac_D[e]'},0,200);
model = addExchangeRxn(model,{'nh4[e]'},-200,0);
model = addExchangeRxn(model,{'o2[e]'},-200,0);
model = addExchangeRxn(model,{'pi[e]'},-200,200);
model = addExchangeRxn(model,{'succ[c]'},0,200);
model = addExchangeRxn(model,{'H2[c]'},0,100);

model=addReaction(model,'glucose_difusion','glc[e] <==> glc[c]');
model=addReaction(model,'GlucosePTS','glc[e] + pep[c] -> g6p[c] + pyr[c]');
model=addReaction(model,'GluKin',{'glc[c]','atp[c]','g6p[c]','adp[c]','h[c]'},[-1 -1 1 1 1],true,0,30);
model=addReaction(model,'PGI','g6p[c] <==> f6p[c]');
model=addReaction(model,'PFK','atp[c] + f6p[c] -> adp[c] + fdp[c] + h[c]');
model=addReaction(model,'ALD','fdp[c] <==> dhap[c] + g3p[c]');
model=addReaction(model,'TPI','dhap[c] <==> g3p[c]');
model=addReaction(model,'GAPDH','g3p[c] + nad[c] + pi[c] <==> 13dpg[c] + h[c] + nadh[c]');
model=addReaction(model,'PGK','13dpg[c] + adp[c] <==> 3pg[c] + atp[c]');
model=addReaction(model,'PGlycerateM','3pg[c] <==> 2pg[c]');
model=addReaction(model,'ENO','2pg[c] <==> h2o[c] + pep[c]');
model=addReaction(model,'PYK','adp[c] + h[c] + pep[c] -> atp[c] + pyr[c]');
model=addReaction(model,'PDH','coa[c] + nad[c] + pyr[c] -> accoa[c] + co2[c] + nadh[c]');
model=addReaction(model,'AKGDH','akg[c] + coa[c] + nad[c] -> co2[c] + nadh[c] + succoa[c]');
model=addReaction(model,'Aconitase','cit[c] <==> icit[c]');
model=addReaction(model,'CS','accoa[c] + h2o[c] + oaa[c] -> cit[c] + coa[c] + h[c] ');
model=addReaction(model,'FUM','fum[c] + h2o[c] <==> mal-L[c]');
model=addReaction(model,'ICDH','icit[c] + nadp[c] -> akg[c] + co2[c] + nadph[c] ');
model=addReaction(model,'MDH','mal-L[c] + nad[c] <==> h[c] + nadh[c] + oaa[c]');
model=addReaction(model,'SUCDH','fadh[c] + succ[c] -> fadh2[c]  + fum[c]');
model=addReaction(model,'SCS','adp[c] + pi[c] + succoa[c] <==> atp[c] + coa[c] + succ[c]');
model=addReaction(model,'DHFRd','dhf[c] + h[c] + nadph[c] <==> nadp[c] + thf[c]');
model=addReaction(model,'PPihydrolase','h2o[c] + ppi[c] -> h[c] + 2 pi[c]');
model=addReaction(model,'PEPC','co2[c] + h2o[c] + pep[c] -> h[c] + oaa[c] + pi[c]');
model=addReaction(model,'G6PDH','g6p[c] + nadp[c] -> 6pgl[c] + h[c] + nadph[c]');
model=addReaction(model,'6PGlactonase','6pgl[c] + h2o[c] -> 6pgc[c] + h[c]');
model=addReaction(model,'GND','6pgc[c] + nadp[c] -> co2[c] + nadph[c] + ru5p-D[c]');
model=addReaction(model,'RibIso','r5p[c] <==> ru5p-D[c]');
model=addReaction(model,'RibEpi','ru5p-D[c] <==> xu5p-D[c]');
model=addReaction(model,'TALA','g3p[c] + s7p[c] <==> e4p[c] + f6p[c] ');
model=addReaction(model,'TKT1','r5p[c] + xu5p-D[c] <==> g3p[c] + s7p[c]');
model=addReaction(model,'TKT2','e4p[c] + xu5p-D[c] <==> f6p[c] + g3p[c]');
model=addReaction(model,'ATPase','adp[c] + 4 h[e] + pi[c] <==> atp[c] + 3 h[c] + h2o[c]');
model=addReaction(model,'ETC_FADH2','fadh2[c] + q8[c] <==> fadh[c] + q8h2[c]');
model=addReaction(model,'r37','h[c] + nadph[c] + trdox[c] -> nadp[c] + trdrd[c]');
model=addReaction(model,'ETOH','accoa[c] + 2 h[c] + 2 nadh[c] <==> coa[c] + etoh[c] + 2 nad[c]');
model=addReaction(model,'ACK','ac[c] + atp[c] <==> actp[c] + adp[c]');
model=addReaction(model,'LDH','h[c] + nadh[c] + pyr[c] <==> lac_D[c] + nad[c]');
model=addReaction(model,'PTA','accoa[c] + pi[c] <==> actp[c] + coa[c]');
model=addReaction(model,'FBPase','fdp[c] + h2o[c] -> f6p[c] + pi[c]');
model=addReaction(model,'PEPCK','atp[c] + oaa[c] -> adp[c] + co2[c] + pep[c]');
model=addReaction(model,'ICL','icit[c] -> glx[c] + succ[c]');
model=addReaction(model,'FUMDHq8','fum[c] + q8h2[c] -> succ[c] + q8[c]');
model=addReaction(model,'PFL','pyr[c] + coa[c] -> accoa[c] + for[c]');
model=addReaction(model,'Hydrogenase','for[c] + h[c] -> co2[c] + H2[c]');
model=addReaction(model,'pntAB','2 h[e] + nadp[c] + nadh[c] -> 2 h[c] + nad[c] + nadph[c]');
model=addReaction(model,'UdhA','nadph[c] + nad[c] -> nadh[c] + nadp[c]');
model=addReaction(model,'ED1',{'6pgc[c]','h2o[c]','KDPG[c]'},[-1 1 1],true,0,20);
model=addReaction(model,'ED2',{'KDPG[c]','pyr[c]','g3p[c]'},[-1 1 1],true,0,20);
model=addReaction(model,'PGlucoseM','g1p[c] <==> g6p[c]');
model=addReaction(model,'GlyK','atp[c] + glyc[c] -> adp[c] + glyc3p[c] + h[c]');
model=addReaction(model,'r278','glyc3p[c] + q8[c] -> dhap[c] + q8h2[c]');
model=addReaction(model,'MALS','accoa[c] + glx[c] + h2o[c] -> coa[c] + h[c] + mal-L[c]');
model=addReaction(model,'ACS','ac[c] + atp[c] + coa[c] -> accoa[c] + amp[c] + ppi[c]');
model=addReaction(model,'r11','atp[c] + g1p[c] + h[c] -> adpglc[c] + ppi[c] ');
model=addReaction(model,'r12','adpglc[c] -> adp[c] + glycogen[c] + h[c] ');
model=addReaction(model,'r42','ru5p-D[c] <==> ara5p[c]');
model=addReaction(model,'r43','2dr5p[c] -> acald[c] + g3p[c]');
model=addReaction(model,'Glyc3PDH','dhap[c] + h[c] + nadh[c] -> glyc3p[c] + nad[c]');
model=addReaction(model,'r46','r1p[c] <==> r5p[c]');
model=addReaction(model,'r47','2dr1p[c] <==> 2dr5p[c]');
model=addReaction(model,'r48','acald[e] <==> acald[c]');
model=addReaction(model,'r49','ac[e] <==> ac[c]');
model=addReaction(model,'r50','nh4[e] <==> nh4[c]');
model=addReaction(model,'r51','co2[e] <==> co2[c]');
model=addReaction(model,'r53','lac_D[e] <==> lac_D[c]');
model=addReaction(model,'r54','etoh[c] <==> etoh[e]');
model=addReaction(model,'r55','for[e] <==> for[c]');
model=addReaction(model,'r56','glyc[c] <==> glyc[e]');
model=addReaction(model,'r57','h2o[e] <==> h2o[c]');
model=addReaction(model,'r58','o2[e] <==> o2[c]');
model=addReaction(model,'r59','pi[e] <==> pi[c]');
model=addReaction(model,'r60','succ[c] -> succ[e]');
model=addReaction(model,'r61','atp[c] + h2o[c] + so4[e] -> adp[c] + h[c] + pi[c] + so4[c]');
model=addReaction(model,'r62','ala-L[c] <==> ala-D[c]');
model=addReaction(model,'r63','asp-L[c] + atp[c] + nh4[c] -> amp[c] + asn-L[c] + h[c] + ppi[c]');
model=addReaction(model,'r64','akg[c] + asp-L[c] <==> glu-L[c] + oaa[c]');
model=addReaction(model,'r65','akg[c] + ala-L[c] <==> glu-L[c] + pyr[c]');
model=addReaction(model,'r66','dkmpp[c] + 3 h2o[c] -> 2kmb[c] + for[c] + 6 h[c] + pi[c]');
model=addReaction(model,'r67','2kmb[c] + glu-L[c] -> akg[c] + met-L[c]');
model=addReaction(model,'r68','5mdru1p[c] -> dkmpp[c] + h2o[c]');
model=addReaction(model,'r69','5mtr[c] + atp[c] -> 5mdr1p[c] + adp[c] + h[c]');
model=addReaction(model,'r70','5mdr1p[c] <==> 5mdru1p[c]');
model=addReaction(model,'r71','acglu[c] + atp[c] -> acg5p[c] + adp[c]');
model=addReaction(model,'r72','acorn[c] + h2o[c] -> ac[c] + orn[c]');
model=addReaction(model,'r73','acorn[c] + akg[c] <==> acg5sa[c] + glu-L[c]');
model=addReaction(model,'r74','amet[c] + h[c] <==> ametam[c] + co2[c]');
model=addReaction(model,'r75','argsuc[c] <==> arg-L[c] + fum[c]');
model=addReaction(model,'r76','asp-L[c] + atp[c] + citr-L[c] -> amp[c] + argsuc[c] + h[c] + ppi[c]');
model=addReaction(model,'r77','2 atp[c] + gln-L[c] + h2o[c] + hco3[c] -> 2 adp[c] + cbp[c] + glu-L[c] + 2 h[c] + pi[c]');
model=addReaction(model,'r78','atp[c] + glu-L[c] -> adp[c] + glu5p[c]');
model=addReaction(model,'r79','glu5p[c] + h[c] + nadph[c] -> glu5sa[c] + nadp[c] + pi[c]');
model=addReaction(model,'r80','glu5sa[c] -> 1pyr5c[c] + h[c] + h2o[c]');
model=addReaction(model,'r81','5mta[c] + h2o[c] -> 5mtr[c] + ade[c]');
model=addReaction(model,'r82','acg5sa[c] + nadp[c] + pi[c] <==> acg5p[c] + h[c] + nadph[c]');
model=addReaction(model,'r83','accoa[c] + glu-L[c] -> acglu[c] + coa[c] + h[c]');
model=addReaction(model,'r84','cbp[c] + orn[c] <==> citr-L[c] + h[c] + pi[c]');
model=addReaction(model,'r85','h[c] + orn[c] -> co2[c] + ptrc[c]');
model=addReaction(model,'r86','1pyr5c[c] + 2 h[c] + nadph[c] -> nadp[c] + pro-L[c]');
model=addReaction(model,'r87','ametam[c] + ptrc[c] -> 5mta[c] + h[c] + spmd[c]');
model=addReaction(model,'r88','h2o[c] + pap[c] -> amp[c] + pi[c]');
model=addReaction(model,'r89','aps[c] + atp[c] -> adp[c] + h[c] + paps[c]');
model=addReaction(model,'r90','acser[c] + h2s[c] -> ac[c] + cys-L[c] + h[c]');
model=addReaction(model,'r91','paps[c] + trdrd[c] -> 2 h[c] + pap[c] + so3[c] + trdox[c]');
model=addReaction(model,'r92','accoa[c] + ser-L[c] <==> acser[c] + coa[c]');
model=addReaction(model,'r93','atp[c] + gtp[c] + h2o[c] + so4[c] -> aps[c] + gdp[c] + pi[c] + ppi[c]');
model=addReaction(model,'r94','3 h2o[c] + h2s[c] + 3 nadp[c] <==> 5 h[c] + 3 nadph[c] + so3[c]');
model=addReaction(model,'r95','glu-L[c] + h2o[c] + nadp[c] <==> akg[c] + h[c] + nadph[c] + nh4[c]');
model=addReaction(model,'r96','atp[c] + glu-L[c] + nh4[c] -> adp[c] + gln-L[c] + h[c] + pi[c]');
model=addReaction(model,'r97','ser-L[c] + thf[c] -> gly[c] + h2o[c] + mlthf[c]');
model=addReaction(model,'r98','3pg[c] + nad[c] -> 3php[c] + h[c] + nadh[c]');
model=addReaction(model,'r99','h2o[c] + pser-L[c] -> pi[c] + ser-L[c]');
model=addReaction(model,'r100','3php[c] + glu-L[c] -> akg[c] + pser-L[c]');
model=addReaction(model,'r101','prfp[c] -> prlp[c]');
model=addReaction(model,'r102','atp[c] + prpp[c] -> ppi[c] + prbatp[c]');
model=addReaction(model,'r103','h2o[c] + histd[c] + 2 nad[c] -> 3 h[c] + his-L[c] + 2 nadh[c]');
model=addReaction(model,'r104','h2o[c] + hisp[c] -> histd[c] + pi[c]');
model=addReaction(model,'r105','glu-L[c] + imacp[c] -> akg[c] + hisp[c]');
model=addReaction(model,'r106','gln-L[c] + prlp[c] -> aicar[c] + eig3p[c] + glu-L[c] + h[c]');
model=addReaction(model,'r107','eig3p[c] -> h2o[c] + imacp[c]');
model=addReaction(model,'r108','h2o[c] + prbamp[c] -> prfp[c]');
model=addReaction(model,'r109','h2o[c] + prbatp[c] -> h[c] + ppi[c] + prbamp[c]');
model=addReaction(model,'r110','atp[c] + r5p[c] <==> amp[c] + h[c] + prpp[c]');
model=addReaction(model,'r111','cyst-L[c] + h2o[c] -> hcys-L[c] + nh4[c] + pyr[c]');
model=addReaction(model,'r112','hom-L[c] + succoa[c] -> coa[c] + suchms[c]');
model=addReaction(model,'r113','atp[c] + h2o[c] + met-L[c] -> amet[c] + pi[c] + ppi[c]');
model=addReaction(model,'r114','5mthf[c] + hcys-L[c] -> met-L[c] + thf[c]');
model=addReaction(model,'r115','cys-L[c] + suchms[c]  -> cyst-L[c] + h[c] + succ[c]');
model=addReaction(model,'r116','asp-L[c] + atp[c] <==> 4pasp[c]  + adp[c]');
model=addReaction(model,'r117','aspsa[c] + nadp[c] + pi[c] <==> 4pasp[c] + h[c] + nadph[c]');
model=addReaction(model,'r118','26dap-M[c] + h[c] -> co2[c] + lys-L[c]');
model=addReaction(model,'r119','26dap-LL[c]  <==> 26dap-M[c]');
model=addReaction(model,'r120','23dhdp[c] + h[c] + nadph[c] -> nadp[c] + thdp[c]');
model=addReaction(model,'r121','aspsa[c] + pyr[c] -> 23dhdp[c]  + h[c] + 2 h2o[c]');
model=addReaction(model,'r122','hom-L[c] + nadp[c] <==> aspsa[c] + h[c] + nadph[c]');
model=addReaction(model,'r123','h2o[c] + sl26da[c] -> 26dap-LL[c] + succ[c]');
model=addReaction(model,'r124','akg[c] + sl26da[c] <==> glu-L[c] + sl2a6o[c]');
model=addReaction(model,'r125','h2o[c] + succoa[c] + thdp[c] -> coa[c] + sl2a6o[c]');
model=addReaction(model,'r126','thr-L[c] <==> acald[c] + gly[c]');
model=addReaction(model,'r127','3dhq[c] <==> 3dhsk[c] + h2o[c]');
model=addReaction(model,'r128','2dda7p[c] -> 3dhq[c] + pi[c]');
model=addReaction(model,'r129','e4p[c] + h2o[c] + pep[c] -> 2dda7p[c] + pi[c]');
model=addReaction(model,'r130','pep[c] + skm5p[c] <==> 3psme[c] + pi[c]');
model=addReaction(model,'r131','anth[c] + prpp[c] -> ppi[c] + pran[c]');
model=addReaction(model,'r132','chor[c] + gln-L[c] -> anth[c] + glu-L[c] + h[c] + pyr[c]');
model=addReaction(model,'r133','chor[c] -> pphn[c]');
model=addReaction(model,'r134','3psme[c] -> chor[c] + pi[c]');
model=addReaction(model,'r135','2cpr5p[c] + h[c] -> 3ig3p[c] + co2[c] + h2o[c]');
model=addReaction(model,'r136','akg[c] + phe-L[c] <==> glu-L[c] + phpyr[c]');
model=addReaction(model,'r137','pran[c] -> 2cpr5p[c]');
model=addReaction(model,'r138','h[c] + pphn[c] -> co2[c] + h2o[c] + phpyr[c]');
model=addReaction(model,'r139','nad[c] + pphn[c] -> 34hpp[c] + co2[c] + nadh[c]');
model=addReaction(model,'r140','3dhsk[c] + h[c] + nadph[c] <==> nadp[c] + skm[c]');
model=addReaction(model,'r141','atp[c] + skm[c] -> adp[c] + h[c] + skm5p[c]');
model=addReaction(model,'r142','3ig3p[c] -> g3p[c] + indole[c]');
model=addReaction(model,'r143','h2o[c] + trp-L[c] <==> indole[c] + nh4[c] + pyr[c]');
model=addReaction(model,'r144','akg[c] + tyr-L[c] <==> 34hpp[c] + glu-L[c]');
model=addReaction(model,'r145','2obut[c] + h[c] + pyr[c] -> 2ahbut[c] + co2[c]');
model=addReaction(model,'r146','2ippm[c] + h2o[c] <==> 3c3hmp[c]');
model=addReaction(model,'r147','3mob[c] + accoa[c] + h2o[c] -> 3c3hmp[c] + coa[c] + h[c]');
model=addReaction(model,'r148','3c4mop[c] + h[c] -> 4mop[c] + co2[c]');
model=addReaction(model,'r149','3c2hmp[c] <==> 2ippm[c] + h2o[c]');
model=addReaction(model,'r150','3c2hmp[c] + nad[c] -> 3c4mop[c] + h[c] + nadh[c]');
model=addReaction(model,'r151','alac-S[c] + h[c] + nadph[c] -> 23dhmb[c] + nadp[c]');
model=addReaction(model,'r152','h[c] + 2 pyr[c] -> alac-S[c] + co2[c]');
model=addReaction(model,'r153','23dhmb[c] -> 3mob[c] + h2o[c]');
model=addReaction(model,'r154','23dhmp[c] -> 3mop[c] + h2o[c]');
model=addReaction(model,'r155','akg[c] + ile-L[c] <==> 3mop[c] + glu-L[c]');
model=addReaction(model,'r156','2ahbut[c] + h[c] + nadph[c] -> 23dhmp[c] + nadp[c]');
model=addReaction(model,'r157','4mop[c] + glu-L[c] -> akg[c] + leu-L[c]');
model=addReaction(model,'r158','thr-L[c] -> 2obut[c] + nh4[c]');
model=addReaction(model,'r159','akg[c] + val-L[c] <==> 3mob[c] + glu-L[c]');
model=addReaction(model,'r160','ara5p[c] + h2o[c] + pep[c] -> kdo8p[c] + pi[c]');
model=addReaction(model,'r161','ckdo[c] + lipidA[c] -> cmp[c] + h[c] + kdolipid4[c]');
model=addReaction(model,'r162','ckdo[c] + kdolipid4[c] -> cmp[c] + h[c] + kdo2lipid4[c]');
model=addReaction(model,'r163','ctp[c] + kdo[c] -> ckdo[c] + ppi[c]');
model=addReaction(model,'r164','h2o[c] + kdo8p[c] -> kdo[c] + pi[c]');
model=addReaction(model,'r165','ACP[c] + atp[c] + ttdcea[c] -> amp[c] + ppi[c] + tdeACP[c]');
model=addReaction(model,'r166','ACP[c] + atp[c] + hdca[c] -> amp[c] + palmACP[c] + ppi[c]');
model=addReaction(model,'r167','ACP[c] + atp[c] + hdcea[c] -> amp[c] + hdeACP[c] + ppi[c]');
model=addReaction(model,'r168','ACP[c] + atp[c] + ocdcea[c] -> amp[c] + octeACP[c] + ppi[c]');
model=addReaction(model,'r169','ACP[c] + atp[c] + ttdca[c] -> amp[c] + myrsACP[c] + ppi[c]');
model=addReaction(model,'r170','adphep-D,D[c] -> adphep-L,D[c]');
model=addReaction(model,'r171','2 ala-D[c] + atp[c] <==> adp[c] + alaala[c] + h[c] + pi[c]');
model=addReaction(model,'r172','atp[c] + gmhep1p[c] + h[c] -> adphep-D,D[c] + ppi[c]');
model=addReaction(model,'r173','gmhep17bp[c] + h2o[c] -> gmhep1p[c] + pi[c]');
model=addReaction(model,'r174','atp[c] + gmhep7p[c] -> adp[c] + gmhep17bp[c] + h[c]');
model=addReaction(model,'r175','0.02 12dgr_EC[c] + atp[c] -> adp[c] + h[c] + 0.02 pa_EC[c]');
model=addReaction(model,'r176','ddcaACP[c] + kdo2lipid4[c] -> ACP[c] + kdo2lipid4L[c]');
model=addReaction(model,'r177','kdo2lipid4L[c] + myrsACP[c] -> ACP[c] + lipa[c]');
model=addReaction(model,'r178','cmp[c] + h[c] + 0.02 pe_EC[c] <==> 0.02 12dgr_EC[c] + cdpea[c]');
model=addReaction(model,'r179','accoa[c] + gam1p[c] -> acgam1p[c] + coa[c] + h[c]');
model=addReaction(model,'r180','glu-D[c] <==> glu-L[c]');
model=addReaction(model,'r181','f6p[c] + gln-L[c] -> gam6p[c] + glu-L[c]');
model=addReaction(model,'r182','g3pe[c] + h2o[c] -> etha[c] + glyc3p[c] + h[c]');
model=addReaction(model,'r183','g3pg[c] + h2o[c] -> glyc[c] + glyc3p[c] + h[c]');
model=addReaction(model,'r184','lipidX[c] + u23ga[c] -> h[c] + lipidAds[c] + udp[c]');
model=addReaction(model,'r185','3 adphep-L,D[c] + 2 cdpea[c] + 3 ckdo[c] + lipa[c] + 2 udpg[c] -> 3 adp[c] + 2 cdp[c] + 3 cmp[c] + 10 h[c] + lps_EC[c] + 2 udp[c]');
model=addReaction(model,'r186','0.02 agpe_EC[c] + h2o[c] -> g3pe[c] + h[c] + 0.36 hdca[c] + 0.07 hdcea[c] + 0.5 ocdcea[c] + 0.02 ttdca[c] + 0.05 ttdcea[c]');
model=addReaction(model,'r187','0.02 agpg_EC[c] + h2o[c] -> g3pg[c] + h[c] + 0.36 hdca[c] + 0.07 hdcea[c] + 0.5 ocdcea[c] + 0.02 ttdca[c] + 0.05 ttdcea[c]');
model=addReaction(model,'r188','uaagmda[c] -> h[c] + peptido_EC[c] + udcpdp[c]');
model=addReaction(model,'r189','h2o[c] + 0.02 pe_EC[c] -> 0.02 agpe_EC[c] + h[c] + 0.36 hdca[c] + 0.07 hdcea[c] + 0.5 ocdcea[c] + 0.02 ttdca[c] + 0.05 ttdcea[c]');
model=addReaction(model,'r190','gam1p[c] <==> gam6p[c]');
model=addReaction(model,'r191','h2o[c] + 0.02 pg_EC[c] -> 0.02 agpg_EC[c] + h[c] + 0.36 hdca[c] + 0.07 hdcea[c] + 0.5 ocdcea[c] + 0.02 ttdca[c] + 0.05 ttdcea[c]');
model=addReaction(model,'r192','udcpp[c] + ugmda[c] -> uagmda[c] + ump[c]');
model=addReaction(model,'r193','s7p[c] -> gmhep7p[c]');
model=addReaction(model,'r194','atp[c] + lipidAds[c] -> adp[c] + h[c] + lipidA[c]');
model=addReaction(model,'r195','3hmrsACP[c] + u3hga[c] -> ACP[c] + h[c] + u23ga[c]');
model=addReaction(model,'r196','h2o[c] + u3aga[c] -> ac[c] + u3hga[c]');
model=addReaction(model,'r197','h[c] + nadph[c] + uaccg[c] -> nadp[c] + uamr[c]');
model=addReaction(model,'r198','pep[c] + uacgam[c] -> pi[c] + uaccg[c]');
model=addReaction(model,'r199','3hmrsACP[c] + uacgam[c] <==> ACP[c] + u3aga[c]');
model=addReaction(model,'r200','acgam1p[c] + h[c] + utp[c] -> ppi[c] + uacgam[c]');
model=addReaction(model,'r201','uacgam[c] + uagmda[c] -> h[c] + uaagmda[c] + udp[c]');
model=addReaction(model,'r202','ala-L[c] + atp[c] + uamr[c] -> adp[c] + h[c] + pi[c] + uama[c]');
model=addReaction(model,'r203','atp[c] + glu-D[c] + uama[c] -> adp[c] + h[c] + pi[c] + uamag[c]');
model=addReaction(model,'r204','26dap-M[c]  + atp[c] + uamag[c] -> adp[c] + h[c] + pi[c] + ugmd[c]');
model=addReaction(model,'r205','alaala[c] + atp[c] + ugmd[c] -> adp[c] + h[c] + pi[c] + ugmda[c]');
model=addReaction(model,'r206','h2o[c] + u23ga[c] -> 2 h[c] + lipidX[c] + ump[c]');
model=addReaction(model,'r207','h2o[c] + udcpdp[c] -> h[c] + pi[c] + udcpp[c]');
model=addReaction(model,'r208','g1p[c] + h[c] + utp[c] <==> ppi[c] + udpg[c]');
model=addReaction(model,'r209','h[c] + mlthf[c] + nadh[c] -> 5mthf[c] + nad[c]');
model=addReaction(model,'r210','10fthf[c] + h2o[c] -> for[c] + h[c] + thf[c]');
model=addReaction(model,'r211','h2o[c] + methf[c] <==> 10fthf[c]');
model=addReaction(model,'r212','mlthf[c] + nadp[c] <==> h[c] + methf[c] + nadph[c]');
model=addReaction(model,'r213','ddcaACP[c] + 2 h[c] + malACP[c] + nadph[c] -> 3hmrsACP[c] + ACP[c] + co2[c] + nadp[c]');
model=addReaction(model,'r214','accoa[c] + atp[c] + hco3[c] <==> adp[c] + h[c] + malcoa[c] + pi[c]');
model=addReaction(model,'r215','accoa[c] + h[c] + malACP[c] -> actACP[c] + co2[c] + coa[c]');
model=addReaction(model,'r216','ctp[c] + h[c] + 0.02 pa_EC[c] <==> 0.02 cdpdag1[c] + ppi[c]');
model=addReaction(model,'r217','actACP[c] + 14 h[c] + 4 malACP[c] + 10 nadph[c] -> 4 ACP[c] + 4 co2[c] + ddcaACP[c] + 5 h2o[c] + 10 nadp[c]');
model=addReaction(model,'r218','actACP[c] + 17 h[c] + 5 malACP[c] + 12 nadph[c] -> 5 ACP[c] + 5 co2[c] + 6 h2o[c] + myrsACP[c] + 12 nadp[c]');
model=addReaction(model,'r219','actACP[c] + 20 h[c] + 6 malACP[c] + 14 nadph[c] -> 6 ACP[c] + 6 co2[c] + 7 h2o[c] + 14 nadp[c] + palmACP[c]');
model=addReaction(model,'r220','actACP[c] + 19 h[c] + 6 malACP[c] + 13 nadph[c] -> 6 ACP[c] + 6 co2[c] + 7 h2o[c] + hdeACP[c] + 13 nadp[c]');
model=addReaction(model,'r221','actACP[c] + 22 h[c] + 7 malACP[c] + 15 nadph[c] -> 7 ACP[c] + 7 co2[c] + 8 h2o[c] + 15 nadp[c] + octeACP[c]');
model=addReaction(model,'r222','ACP[c] + malcoa[c] <==> coa[c] + malACP[c]');
model=addReaction(model,'r223','glyc3p[c] + 0.14 hdeACP[c] + 0.04 myrsACP[c] + octeACP[c] + 0.72 palmACP[c] + 0.1 tdeACP[c] -> 2 ACP[c] + 0.02 pa_EC[c]');
model=addReaction(model,'r224','h2o[c] + 0.02 pgp_EC[c] -> 0.02 pg_EC[c] + pi[c]');
model=addReaction(model,'r225','0.02 cdpdag1[c] + glyc3p[c] <==> cmp[c] + h[c] + 0.02 pgp_EC[c]');
model=addReaction(model,'r226','h[c] + 0.02 ps_EC[c] -> co2[c] + 0.02 pe_EC[c]');
model=addReaction(model,'r227','0.02 cdpdag1[c] + ser-L[c] <==> cmp[c] + h[c] + 0.02 ps_EC[c]');
model=addReaction(model,'r228','damp[c] + h2o[c] -> dad-2[c] + pi[c]');
model=addReaction(model,'r229','adn[c] + atp[c] -> adp[c] + amp[c] + h[c]');
model=addReaction(model,'r230','amp[c] + atp[c] <==> 2 adp[c]');
model=addReaction(model,'r231','atp[c] + cmp[c] <==> adp[c] + cdp[c]');
model=addReaction(model,'r232','atp[c] + dcmp[c] <==> adp[c] + dcdp[c]');
model=addReaction(model,'r233','atp[c] + damp[c] <==> adp[c] + dadp[c]');
model=addReaction(model,'r234','atp[c] + dgmp[c] <==> adp[c] + dgdp[c]');
model=addReaction(model,'r235','atp[c] + gmp[c] <==> adp[c] + gdp[c]');
model=addReaction(model,'r236','atp[c] + cdp[c] <==> adp[c] + ctp[c]');
model=addReaction(model,'r237','atp[c] + dudp[c] <==> adp[c] + dutp[c]');
model=addReaction(model,'r238','atp[c] + gdp[c] <==> adp[c] + gtp[c]');
model=addReaction(model,'r239','atp[c] + udp[c] <==> adp[c] + utp[c]');
model=addReaction(model,'r240','adn[c] + pi[c] <==> ade[c] + r1p[c]');
model=addReaction(model,'r241','dad-2[c] + pi[c] <==> 2dr1p[c] + ade[c]');
model=addReaction(model,'r242','adp[c] + trdrd[c] -> dadp[c] + h2o[c] + trdox[c]');
model=addReaction(model,'r243','cdp[c] + trdrd[c] -> dcdp[c] + h2o[c] + trdox[c]');
model=addReaction(model,'r244','gdp[c] + trdrd[c] -> dgdp[c] + h2o[c] + trdox[c]');
model=addReaction(model,'r245','trdrd[c] + utp[c] -> dutp[c] + h2o[c] + trdox[c]');
model=addReaction(model,'r246','dump[c] + mlthf[c] -> dhf[c] + dtmp[c]');
model=addReaction(model,'r247','atp[c] + ump[c] <==> adp[c] + udp[c]');
model=addReaction(model,'r248','atp[c] + dump[c] <==> adp[c] + dudp[c]');
model=addReaction(model,'r249','25aics[c] <==> aicar[c] + fum[c]');
model=addReaction(model,'r250','asp-L[c] + gtp[c] + imp[c] -> dcamp[c] + gdp[c] + 2 h[c] + pi[c]');
model=addReaction(model,'r251','dcamp[c] <==> amp[c] + fum[c]');
model=addReaction(model,'r252','asp-L[c] + cbp[c] -> cbasp[c] + h[c] + pi[c]');
model=addReaction(model,'r253','atp[c] + gln-L[c] + h2o[c] + utp[c] -> adp[c] + ctp[c] + glu-L[c] + 2 h[c] + pi[c]');
model=addReaction(model,'r254','dhor-S[c] + q8[c] -> orot[c] + q8h2[c]');
model=addReaction(model,'r255','dhor-S[c] + h2o[c] <==> cbasp[c] + h[c]');
model=addReaction(model,'r256','gln-L[c] + h2o[c] + prpp[c] -> glu-L[c] + ppi[c] + pram[c]');
model=addReaction(model,'r257','atp[c] + gln-L[c] + h2o[c] + xmp[c] -> amp[c] + glu-L[c] + gmp[c] + 2 h[c] + ppi[c]');
model=addReaction(model,'r258','h2o[c] + imp[c] <==> fprica[c]');
model=addReaction(model,'r259','h2o[c] + imp[c] + nad[c] -> h[c] + nadh[c] + xmp[c]');
model=addReaction(model,'r260','orot5p[c] + ppi[c] <==> orot[c] + prpp[c]');
model=addReaction(model,'r261','h[c] + orot5p[c] -> co2[c] + ump[c]');
model=addReaction(model,'r262','air[c] + atp[c] + hco3[c] -> 5caiz[c] + adp[c] + h[c] + pi[c]');
model=addReaction(model,'r263','5aizc[c] <==> 5caiz[c]');
model=addReaction(model,'r264','atp[c] + fpram[c] -> adp[c] + air[c] + 2 h[c] + pi[c]');
model=addReaction(model,'r265','10fthf[c] + aicar[c] <==> fprica[c] + thf[c]');
model=addReaction(model,'r266','5aizc[c] + asp-L[c] + atp[c] <==> 25aics[c] + adp[c] + h[c] + pi[c]');
model=addReaction(model,'r267','atp[c] + fgam[c] + gln-L[c] + h2o[c] -> adp[c] + fpram[c] + glu-L[c] + h[c] + pi[c]');
model=addReaction(model,'r268','10fthf[c] + gar[c] <==> fgam[c] + h[c] + thf[c]');
model=addReaction(model,'r269','atp[c] + gly[c] + pram[c] <==> adp[c] + gar[c] + h[c] + pi[c]');
model=addReaction(model,'r271','co2[c] + h2o[c] <==> h[c] + hco3[c]');
model=addReaction(model,'Protein_amino_acids','0.113 ala-L[c] + 0.0512 arg-L[c] + 0.0532 asn-L[c] + 0.0532 asp-L[c] + 0.0176 cys-L[c] + 0.0599 gln-L[c] + 0.0599 glu-L[c] + 0.0872 gly[c] + 0.0182 his-L[c] + 0.0493 ile-L[c] + 0.0541 leu-L[c] + 0.0605 lys-L[c] + 0.0259 met-L[c] + 0.0350 phe-L[c] + 0.0416 pro-L[c] + 0.0501 ser-L[c] + 0.0545 thr-L[c] + 0.0114 trp-L[c] + 0.0290 tyr-L[c] + 0.0752 val-L[c] -> 1 Aaprot[c]');
model=addReaction(model,'Protein_translation','1 Aaprot[c] + 4 atp[c] + 3 h2o[c] -> 4 adp[c] + 4.77 Biom_Prot[c] + 4 h[c] + 4 pi[c]');
model=addReaction(model,'DNA','2 atp[c] + 0.246 damp[c] + 0.254 dcmp[c] + 0.254 dgmp[c] + 0.246 dtmp[c] + h2o[c] -> 2 adp[c] + 9.75 DNA[c] + 2 h[c] + 2 pi[c]');
model=addReaction(model,'RNA','0.262 amp[c] + 2 atp[c] + 0.2 cmp[c] + 0.322 gmp[c] + h2o[c] + 0.216 ump[c] -> 2 adp[c] + 2 h[c] + 2 pi[c] + 9.58 RNAtot[c]');
model=addReaction(model,'Growth','0.7400 Biom_Prot[c] + 0.00865 DNA[c] + 0.0109 etha[c] + 0.0122 glyc[c] + 0.000956 glycogen[c] + 0.00304 hdca[c] + 0.00235 hdcea[c] + 0.0000543 lps_EC[c] + 0.00154 ocdcea[c] + 0.000156 peptido_EC[c] + 0.000784 ptrc[c] + 0.0535 RNAtot[c] + 0.000159 spmd[c] -> 1 BuildingBlocks');

% Energetic parameters 
Kx=0.46;%+/-0.27
mATP=0.075;%+/-0.015
d=1.49;%+/-0.26

model=addReaction(model,'ATPM','atp[c] + h2o[c] -> adp[c] + h[c] + pi[c]');
model=changeRxnBounds(model,'ATPM',mATP,'l');
model=addReaction(model,'Growth_coupled_ATP',{'BuildingBlocks','atp[c]','h2o[c]','Biomass','adp[c]','pi[c]','h[c]'},[-1 -Kx -Kx 1 Kx Kx Kx],true,0,100);
model=addReaction(model,'ETC_O2',{'h[c]','o2[c]','q8h2[c]','h[e]','h2o[c]','q8[c]'},[-2*d -0.5 -1 2*d 1 1],true,0,20);
model=addReaction(model,'ETC_NADH',{'h[c]','nadh[c]','h[c]','q8[c]','h[e]','nad[c]','q8h2[c]'},[-2*d -1 -1 -1 2*d 1 1],true,0,20);
model=addReaction(model,'ETC_FOR',{'for[e]','h[c]','q8[c]','co2[c]','h[e]','q8h2[c]'},[-1 -1*(d+1) -1 1 d 1],true,0,20);

model = addExchangeRxn(model,{'Biomass'},0,1000);

% conversion factor from Cmol to gDCW
Mw=24.518;
%VERY IMPORTANT! The constrainst should be introduced in units of mol/Cmol_of_Biomass/h.
% For example, q_glucose of 14 mmol/gCDW/h is represented as -14*(Mw/1000)
model = changeRxnBounds(model,'EX_glc[e]',-8*(Mw/1000),'l');
model = changeRxnBounds(model,'EX_o2[e]',-15*(Mw/1000),'l');
%To enable anaerobic growth
model = changeRxnBounds(model,'ATPase',-20*(Mw/1000),'l');

model = changeObjective(model,'EX_Biomass');
% FBAsolution = optimizeCbModel(model,'max');
% printFluxVector(model, FBAsolution.x,1,1);

clear d Kx mATP ans FBAsolution

save cellular_model
model_sbml = writeCbModel(model,'format','sbml');