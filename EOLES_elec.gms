$Ontext
This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
Written by Behrang Shirizadeh, October 2020
$Offtext
*-------------------------------------------------------------------------------
*                                Defining the sets
*-------------------------------------------------------------------------------
sets     h                                               /0*8759/
         first(h)        'first hour'
         last(h)         'last hour'
         m               'month'                         /1*12/
         tec             'technology'                    /offshore, onshore, pv, river, lake, biogas1, biogas2, ocgt, ccgt-ccs, ngas1, ngas2, nuc, phs, battery, methanation1, methanation2/
         gen(tec)        'power plants'                  /offshore, onshore, pv, river, lake, biogas1, biogas2, ocgt, ccgt-ccs, ngas1, ngas2, nuc/
         vre(tec)        'variable tecs'                 /offshore, onshore, pv/
         ncomb(tec)      'non-combustible generation'    /offshore, onshore, pv, river, lake, nuc, phs, battery/
         str(tec)        'storage technologies'          /phs, battery, methanation1, methanation2/
         frr(tec)        'technologies for upward FRR'   /lake, phs, battery, ocgt, ccgt-ccs, nuc/
         scenEPR         'EPR cost scenarios'            /1*7/
         scenCO2         'CO2 tax scenarios'             /1*6/
         scenWind        'wind power cost scenarios'     /1*3/
         scenPV          'solar power cost scenarios'    /1*3/
         costtype        'fixed costs cost items'        /capital,operational/
         wind            'wind power technologies'       /offshore,onshore/
         costCO2         'CO2 cost items'                /ngas1,ngas2,biogas2,methanation2/
;
first(h) = ord(h)=1;
last(h) = ord(h)=card(h);
alias(h,hh);
*-------------------------------------------------------------------------------
*                                Inputs
*-------------------------------------------------------------------------------
$Offlisting
parameter month(h) 'relating the hours to the months'
/
$ondelim
$include  inputs/month.csv
$offdelim
/;
parameter load_factor(vre,h) 'Production profiles of VRE'
/
$ondelim
$include  inputs/vre_profiles2006.csv
$offdelim
/;
parameter demand(h) 'demand profile in each hour in GW'
/
$ondelim
$include inputs/demand2050_ademe.csv
$offdelim
/;
Parameter lake_inflows(m) 'monthly lake inflows in GWh'
/
$ondelim
$include  inputs/lake_inflows.csv
$offdelim
/ ;
parameter gene_river(h) 'hourly run of river power generation in GWh'
/
$ondelim
$include  inputs/run_of_river.csv
$offdelim
/ ;
parameter epsilon(vre) 'additional FRR requirement for variable renewable energies because of forecast errors'
/
$ondelim
$include  inputs/reserve_requirements.csv
$offdelim
/ ;
parameter capa_ex(tec) 'existing capacities of the technologies by December 2017 in GW'
/
$ondelim
$include  inputs/existing_capas.csv
$offdelim
/ ;
parameter capa_max(vre) 'maximum capacities of the technologies in GW'
/
$ondelim
$include  inputs/max_capas.csv
$offdelim
/ ;
parameter capex(tec) 'annualized power capex cost in M€/GW/year'
/
$ondelim
$include  inputs/annuitiesccs.csv
$offdelim
/ ;
parameter capex_en(str) 'annualized energy capex cost of storage technologies in M€/GWh/year'
/
$ondelim
$include  inputs/str_annuitiesccs.csv
$offdelim
/ ;
parameter fOM(tec) 'annualized fixed operation and maintenance costs M€/GW/year'
/
$ondelim
$include  inputs/fO&Mccs.csv
$offdelim
/ ;
Parameter vOM(tec) 'Variable operation and maintenance costs in M€/GWh'
/
$ondelim
$include  inputs/vO&Mccs.csv
$offdelim
/ ;
parameter costs1(scenEPR,costtype) 'the values of changing parameters'
/
$ondelim
$include inputs/scenariosEPR.csv
$offdelim
/ ;
parameter costs2(scenCO2,costCO2) 'the values of changing parameters'
/
$ondelim
$include inputs/scenariosCO2.csv
$offdelim
/ ;
parameter costs3(scenWind,wind,costtype) 'the values of gas price in each scenario'
/
$ondelim
$include inputs/scenariosWind.csv
$offdelim
/;
parameter costs4(scenPV,costtype) 'the values of gas price in each scenario'
/
$ondelim
$include inputs/scenariosPV.csv
$offdelim
/;
parameter fixed_costs(tec) 'yearly fixed cost of each tec in M€/GW/year' ;
fixed_costs(tec) = capex(tec) + fOM(tec);
parameter s_capex(str) 'charging related annuity of storage in M€/GW/year' /PHS 0, battery 0, methanation1 84.16086, methanation2 84.16086/;
parameter s_opex(str)    'charging related fOM of storage in M€/GW/year'   /PHS 0, battery 0, methanation1 59.25, methanation2 59.25/;
parameter eta_in(str) 'charging efifciency of storage technologies' /PHS 0.95, battery 0.9, methanation1 0.59, methanation2 0.59/;
parameter eta_out(str) 'discharging efficiency of storage technolgoies' /PHS 0.9, battery 0.95, methanation1 0.45, methanation2 0.45/;
scalar eta_ocgt 'efficiency of OCGT power plants' /0.45/;
scalar eta_ccgt 'efifciency of CCGT power plants with CCS' /0.55/;
scalar cf_nuc 'maximum capacity factor of nuclear power plants' /0.90/;
scalar ramp_rate 'maximum ramp up/down rate for nuclear power plant' /0.25/;
scalar cf_ccgt 'maximum capaity factor of CCGT plant for a year' /0.85/;
scalar pump_capa 'pumping capacity in GW' /9.3/;
scalar max_phs 'maximum volume of energy can be stored in PHS reservoir in TWh' /0.18/;
scalar max_biogas 'maxium energy can be generated by biogas in TWh' /15/;
scalar load_uncertainty 'uncertainty coefficient for hourly demand' /0.01/;
scalar delta 'load variation factor'     /0.1/;
parameter CO2_tax(scenCO2) 'CO2 tax for each tax scenario' /1 0, 2 100, 3 200, 4 300, 5 400, 6 500/;
parameter vOM1(*) 'O&M cost on unit CO2 produced';
parameter vOM2(*) 'transport of CO2 cost for BECCS';
parameter vOM3(*) 'the remuneration of negative emissions from BECCS';
*-------------------------------------------------------------------------------
*                                Model
*-------------------------------------------------------------------------------
variables        GENE(tec,h)     'hourly energy generation in TWh'
                 CAPA(tec)       'overal yearly installed capacity in GW'
                 STORAGE(str,h)  'hourly electricity input of battery storage GW'
                 S(str)          'charging power capacity of each storage technology'
                 STORED(str,h)   'energy stored in each storage technology in GWh'
                 CAPACITY(str)   'energy volume of storage technologies in GWh'
                 RSV(frr,h)      'required upward frequency restoration reserve in GW'
                 E_DAC(h)        'Needed electricity for direct air capturing as GWh/MtCO2'
                 G_DAC(h)        'Needed gas for direct air capturing as GWh/MtCO2'
                 COST            'final investment cost in b€'

positive variables GENE(tec,h),CAPA(tec),STORAGE(str,h), S(str),STORED(str,h),CAPACITY(str),RSV(frr,h),E_DAC(h),G_DAC(h);

equations        gene_vre        'variables renewable profiles generation'
                 gene_capa       'capacity and genration relation for technologies'
                 combustion1     'the relationship of combustible technologies'
                 combustion2     'the relationship of combustible technologies'
                 capa_frr        'capacity needed for the secondary reserve requirements'
                 storing         'the definition of stored energy in the storage options'
                 storage_const   'storage in the first hour is equal to the storage in the last hour'
                 lake_res        'constraint on water for lake reservoirs'
                 stored_cap      'maximum energy that is stored in storage units'
                 storage_capa1   'the capacity with hourly charging relationship of storage'
                 storage_capa2   'storage power limit'
                 biogas_const    'maximum energy can be produced by biogas'
                 nuc_cf          'the yearly capacity factor of nuclear power plants should not pass 80%'
                 nuc_up          'Nuclear power plant upward flexibility flexibility'
                 nuc_down        'Nuclear power plant downward flexibility flexibility'
                 ccgt_cf         'the yearly capacity factor of CCGT'
                 reserves        'FRR requirement'
                 adequacy        'supply/demand relation'
                 obj             'the final objective function which is COST';

gene_vre(vre,h)..                GENE(vre,h)             =e=     CAPA(vre)*load_factor(vre,h);
gene_capa(tec,h)..               CAPA(tec)               =g=     GENE(tec,h);
combustion1(h)..                 GENE('ocgt',h)          =e=     GENE('methanation1',h) + GENE('biogas1',h) + GENE('ngas1',h);
combustion2(h)..                 GENE('ccgt-ccs',h)      =e=     GENE('methanation2',h) + GENE('biogas2',h) + GENE('ngas2',h);
capa_frr(frr,h)..                CAPA(frr)               =g=     GENE(frr,h) + RSV(frr,h);
storing(h,h+1,str)..             STORED(str,h+1)         =e=     STORED(str,h) + STORAGE(str,h)*eta_in(str) - GENE(str,h)/eta_out(str);
storage_const(str,first,last)..  STORED(str,first)       =e=     STORED(str,last) + STORAGE(str,last)*eta_in(str) - GENE(str,last)/eta_out(str);
lake_res(m)..                    lake_inflows(m)         =g=     sum(h$(month(h) = ord(m)),GENE('lake',h))/1000;
stored_cap(str,h)..              STORED(str,h)           =l=     CAPACITY(str);
storage_capa1(str,h)..           S(str)                  =g=     STORAGE(str,h);
storage_capa2(str)..             S(str)                  =l=     CAPA(str);
biogas_const..                   sum(h,GENE('biogas1',h))+sum(h,GENE('biogas2',h)) =l=     max_biogas*1000;
nuc_cf..                         sum(h,GENE('nuc',h))    =l=     CAPA('nuc')*cf_nuc*8760;
nuc_up(h,h+1)..                  GENE('nuc',h+1) + RSV('nuc',h+1) =l= GENE('nuc',h) + ramp_rate*(CAPA('nuc')-GENE('nuc',h))   ;
nuc_down(h,h+1)..                GENE('nuc',h+1) =g= GENE('nuc',h)*(1 - ramp_rate)   ;
ccgt_cf..                        sum(h,GENE('ccgt-ccs',h))=l=    CAPA('ccgt-ccs')*cf_ccgt*8760;
reserves(h)..                    sum(frr, RSV(frr,h))    =e=     sum(vre,epsilon(vre)*CAPA(vre))+ demand(h)*load_uncertainty*(1+delta);
adequacy(h)..                    sum(ncomb,GENE(ncomb,h))+GENE('ocgt',h)+GENE('ccgt-ccs',h)    =g=     demand(h) + sum(str,STORAGE(str,h));
obj..                            COST                    =e=     (sum(tec,(CAPA(tec)-capa_ex(tec))*capex(tec))+ sum(str,CAPACITY(str)*capex_en(str))+sum(tec,(CAPA(tec)*fOM(tec)))+ sum(str,S(str)*(s_capex(str)+s_opex(str))) + sum((tec,h),GENE(tec,h)*vOM(tec)))/1000;
*-------------------------------------------------------------------------------
*                                Initial and fixed values
*-------------------------------------------------------------------------------
CAPA.fx('phs') = pump_capa;
CAPA.fx('river')= capa_ex('river');
CAPA.fx('lake') = 12.855;
GENE.up('river',h) = gene_river(h)*capa_ex('river');
S.fx('phs') = pump_capa;
CAPACITY.fx('phs') = max_phs*1000;
CAPA.up('offshore') = capa_max('offshore');
CAPA.up('onshore') = capa_max('onshore');
CAPA.up('pv') = capa_max('pv');
*-------------------------------------------------------------------------------
*                                Model options
*-------------------------------------------------------------------------------
model EOLES_elec /all/;
*-------------------------------------------------------------------------------
option solvelink=0;
option RESLIM = 1000000;
option lp=CPLEX;
option Savepoint=1;
option solveopt = replace;
option limcol = 0;
option limrow = 0;
option SOLPRINT = OFF;
option solvelink=0;
$onecho > cplex.opt
objrng all
rhsrng all
rngrestart rng.inc
$offecho
EOLES_elec.optfile=1; EOLES_elec.dictfile=2;
*-------------------------------------------------------------------------------
*                                Solve statement
*-------------------------------------------------------------------------------
$If exist EOLES_elec_p.gdx execute_loadpoint 'EOLES_elec_p';
parameter sumdemand      'the whole demand per year in TWh';
parameter gene_tec(tec) 'Overall yearly energy generated by the technology in TWh';
parameter sumgene        'the whole generation per year in TWh';
parameter sum_FRR 'the whole yearly energy budgeted for reserves in TWh';
parameter reserve(frr) 'capacity allocated for reserve from each FRR tech in GW';
parameter lcoe_sys1;
parameter lcoe_sys2;
parameter str_loss 'yearly storage related loss in % of power production';
parameter lc 'load curtailment of the network';
parameter spot_price(h) 'marginal cost'    ;
parameter marginal_cost 'average value over the year of spot price in €/MWh';
parameter CO2_positive 'positive CO2 emission in MtCO2/year';
parameter CO2_negative 'negative CO2 emission in MtCO2/year';
parameter CO2_emission 'the overall CO2 balance in MtCO2/year';
parameter negative_CCS 'yearly CO2 captured by CCS in MtCO2/year';
parameter positive_CCS 'yearly CO2 emitted by CCS in MtCO2/year';
parameter real_cost 'the overall real cost of the system without considering carbon tax or remunerations in b€';
parameter nSTORAGE(str,h);
file summary_csv /'outputs/results.csv' / ;
put summary_csv;
summary_csv.pc=5;
put 'scenCO2','scenEPR','scenWind','scenPV','cost','real_cost','emission','CCS_emission','CCS_capture', 'offshore', 'onshore', 'PV', 'biogas1','biogas2','nuc','ngas1','ngas2', 'battery', 'methanation1', 'methanation2', 'battery_v', 'methanation1_v', 'methanation2_v', 'offshore_gene', 'onshore_gene','PV_gene','biogas1_gene','biogas_2_gene','nuc_gene','ngas1_gene','ngas2_gene', 'battery_gene', 'PHS_gene','methanation1_gene','methanation1_gene','LCOE2', 'marginal_cost', 'LC'/ ;
loop((scenEPR,scenCO2,scenWind,scenPV),
capex('nuc') = costs1(scenEPR,'capital');
fOM('nuc') = costs1(scenEPR,'operational');
vOM('ngas2') = (costs2(scenCO2,'ngas2')+50)/1000;
vOM('ngas1') = (costs2(scenCO2,'ngas1')+50*55/45)/1000;
vOM('biogas2') = costs2(scenCO2,'biogas2')/1000;
vOM('methanation2') = costs2(scenCO2,'methanation2')/1000;
capex('offshore') = costs3(scenWind,'offshore','capital');
capex('onshore') = costs3(scenWind,'onshore','capital');
fOM('offshore') = costs3(scenWind,'offshore','operational');
fOM('onshore') = costs3(scenWind,'onshore','operational');
capex('pv') = costs4(scenPV,'capital');
fOM('pv') = costs4(scenPV,'operational');
Solve EOLES_elec using lp minimizing COST;
*-------------------------------------------------------------------------------
*                                Display statement
*-------------------------------------------------------------------------------
sumdemand =  sum(h,demand(h))/1000;
gene_tec(tec) = sum(h,GENE.l(tec,h))/1000;
sumgene = sum((gen,h),GENE.l(gen,h))/1000 - gene_tec('ocgt')-gene_tec('ccgt-ccs');
sum_FRR = sum((h,frr),RSV.l(frr,h))/1000;
reserve(frr) = smax(h,RSV.l(frr,h));
lcoe_sys1 = cost.l*1000/sumgene;
lcoe_sys2 = cost.l*1000/sumdemand;
str_loss = (sum((str,h),STORAGE.l(str,h))-sum(str,gene_tec(str)*1000))/(sumgene*10);
lc = ((sumgene - sumdemand)*100/sumgene) - str_loss;
spot_price(h) = 1000000*adequacy.m(h);
marginal_cost = sum(h,spot_price(h))/8760;
CO2_positive = ((sum(h,GENE.l('ngas2',h))*40)+(sum(h,GENE.l('ngas1',h))*320))/1000000;
CO2_negative = sum(h,(GENE.l('biogas2',h)+GENE.l('methanation2',h)))*320/1000000;
CO2_emission = CO2_positive - CO2_negative;
negative_CCS =  sum(h,(GENE.l('biogas2',h)+GENE.l('methanation2',h)))*320/1000000  ;
positive_CCS =  CO2_positive;
real_cost = COST.l - (CO2_positive - CO2_negative)*CO2_tax(scenCO2)/1000;
nSTORAGE(str,h) = 0 - STORAGE.l(str,h);
*-------------------------------------------------------------------------------
display cost.l;
display capa.l;
display gene_tec;
display sumdemand; display sumgene;
display lcoe_sys1; display lcoe_sys2;
display CO2_positive; display CO2_negative; display CO2_emission;
display negative_CCS; display positive_CCS;
display CAPACITY.l;
display lc; display str_loss; display marginal_cost;
display real_cost;
*-------------------------------------------------------------------------------
*                                Output
*-------------------------------------------------------------------------------
put summary_csv;
put scenCO2.tl, scenEPR.tl, scenWind.tl, scenPV.tl, COST.l,real_cost,CO2_emission,positive_CCS,negative_CCS,CAPA.l('offshore'),CAPA.l('onshore'),CAPA.l('pv'), CAPA.l('biogas1'),CAPA.l('biogas2'),CAPA.l('nuc'),CAPA.l('ngas1'),CAPA.l('ngas2'),CAPA.l('battery'),CAPA.l('methanation1'), CAPA.l('methanation2'), CAPACITY.l('battery'), CAPACITY.l('methanation1'), CAPACITY.l('methanation2'), gene_tec('offshore'),gene_tec('onshore'),gene_tec('pv'), gene_tec('biogas1'),gene_tec('biogas2'),gene_tec('nuc'),gene_tec('ngas1'),gene_tec('ngas2'),gene_tec('battery'),gene_tec('phs'),gene_tec('methanation1'), gene_tec('methanation2'),lcoe_sys2, marginal_cost, LC/ ;
);
*-------------------------------------------------------------------------------
*                                The End :D
*-------------------------------------------------------------------------------
