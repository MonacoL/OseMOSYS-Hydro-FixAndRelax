# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 12:13:09 2017

@author:     
"""

from __future__ import division
from decimal import *
import pyomo
from pyomo.environ import *
from pyomo.core import *
from pyomo.opt import SolverFactory
import datetime
import numpy as np
import csv
import os, sys
sys.path.insert(0, "../")
try:
    from __init__ import *
except ImportError:
    print('No Import')

import modeling_output.Data as DATA
import automation.telegramsend as TELEGRAM

path = os.path.abspath(sys.argv[0])[:-15]

start=datetime.datetime.now()

#####IMPOSTAZIONI MODELLO#########
PastOrFuture="Future"
pathDatiRun="20190316_1403_2014_2019_fer_DatiRun"
scartoupper=5
scartolower=5
############################## 

#tiro fuori range anni e nome del file di input 
anni = open(path[:-7]+"dat_files/FixAndRelax/"+PastOrFuture+"/"+pathDatiRun+".dat", "r")
vettdatirun = []
for val in anni.readlines():
    vettdatirun.append(val)
anni.close()


VecchieOrNuove=str(vettdatirun[1][:-1])
MinYearToDiscount=int(vettdatirun[3][:-1])
MaxYearToDiscount=int(vettdatirun[4][:-1])
MinGlobalYear=int(vettdatirun[5][:-1])
MaxGlobalYear=int(vettdatirun[6][:-1])
MinYear=int(vettdatirun[7][:-1])
MaxYear=int(vettdatirun[8][:-1])
OverlappingSize=int(vettdatirun[9][:-1])
EmPen=float(vettdatirun[10][:-1])
nomefile=str(vettdatirun[11][:-1])

#########################

model = AbstractModel()
###############
#    Sets     #
############### 

model.YEAR=Set()
model.YEARSMINUSLAST=Set()
model.ACTIVEYEARS=Set()
model.MILPYEARS=Set()
model.ACTIVEMILPYEARS=Set()
if MinGlobalYear<MinYear:
    model.FIXEDMILPYEARS=Set()
model.TECHNOLOGY = Set()
model.TIMESLICE = Set()
model.FUEL = Set()
model.MODE_OF_OPERATION = Set()
model.REGION = Set()
model.SEASON = Set()
model.EMISSION = Set()

if VecchieOrNuove=="Vecchie":
    Reservoirs=DATA.ReservoirsPast()
else:
    Reservoirs=DATA.ReservoirsFuture() 

model.RESERVOIR = Set()
model.OLDRESERVOIR = Set()
model.I = RangeSet(1,len(Reservoirs[-1][11])-1)
#    model.I2 = RangeSet(0,len(Reservoirs[-1][11])-1)
model.J = RangeSet(1,len(Reservoirs[-1][11])-1)


#####################
#    Parameters     #
#####################


########            Global                      #############

model.YearSplit = Param(model.TIMESLICE, model.YEAR)
model.DiscountRate = Param(model.REGION, default=0.05)
model.DaySplit = Param(default=0.00264)
model.Conversionls = Param(model.TIMESLICE, model.SEASON, default=0)
model.DepreciationMethod = Param(model.REGION, default=1)

########            Demands                     #############

model.SpecifiedAnnualDemand = Param(model.REGION, model.FUEL, model.ACTIVEYEARS, default=0)
model.SpecifiedDemandProfile = Param(model.REGION, model.FUEL, model.TIMESLICE, model.ACTIVEYEARS, default=0)

#########           Performance                 #############

model.CapacityToActivityUnit = Param(model.REGION, model.TECHNOLOGY, default=1)
model.CapacityFactor = Param(model.REGION, model.TECHNOLOGY, model.TIMESLICE, model.ACTIVEYEARS, default=1)
model.AvailabilityFactor = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=1)
model.OperationalLife = Param(model.REGION, model.TECHNOLOGY, default=1)
model.ResidualCapacity = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=0)
model.InputActivityRatio = Param(model.REGION, model.TECHNOLOGY, model.FUEL, model.MODE_OF_OPERATION, model.ACTIVEYEARS, default=0)
model.OutputActivityRatio = Param(model.REGION, model.TECHNOLOGY, model.FUEL, model.MODE_OF_OPERATION, model.ACTIVEYEARS, default=0)

#########           Technology Costs            #############

model.CapitalCost = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=0)
model.VariableCost = Param(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.YEAR, default=0.00001)
model.FixedCost = Param(model.REGION, model.TECHNOLOGY, model.YEAR, default=0)

#########           Capacity Constraints        #############

model.CapacityOfOneTechnologyUnit = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=0.0)
model.TotalAnnualMaxCapacity = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=99999)
model.TotalAnnualMinCapacity = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=0)

#########           Investment Constraints      #############

model.TotalAnnualMaxCapacityInvestment = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=99999)
model.TotalAnnualMinCapacityInvestment = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=0)


#########           Activity Constraints        #############

model.TotalTechnologyAnnualActivityUpperLimit = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=99999)
model.TotalTechnologyAnnualActivityLowerLimit = Param(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, default=0)
model.TotalTechnologyModelPeriodActivityUpperLimit = Param(model.REGION, model.TECHNOLOGY, default=99999)
model.TotalTechnologyModelPeriodActivityLowerLimit = Param(model.REGION, model.TECHNOLOGY, default=0)


#########           Hydro Related Parameters            #############
model.RhoWater = Param(default=1000)
model.Gravity = Param(default=9.81)
model.FlowUnittoYearConversion = Param(default=31536000)
model.WattsToModelUnits = Param(default=1e-9)
model.HydroTurbineEfficiency = Param(model.TECHNOLOGY, default=0)
model.TechnologyFromHydro = Param(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.RESERVOIR, default=0)
model.ReservoirHead = Param(model.RESERVOIR,default=0)

model.ReservoirStorageVolume = Param (model.RESERVOIR, default=0)
model.MinimumOperatingLevel = Param (model.RESERVOIR, default=0)
model.ReservoirLevelStart = Param(model.REGION, model.RESERVOIR,default=0)
model.StartRelease = Param(model.RESERVOIR, model.SEASON, default=0)
model.ReservoirExternalInflow = Param(model.REGION, model.SEASON, model.ACTIVEMILPYEARS, model.RESERVOIR,default=0)
model.Evap = Param(model.RESERVOIR, model.SEASON, default=0)
model.DownstreamReservoirTag = Param (model.RESERVOIR, model.RESERVOIR, default=0)
model.InflowLagTime = Param(model.RESERVOIR, model.RESERVOIR, default=1)

model.RuleCurve = Param(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, default=0)
model.MEF = Param(model.RESERVOIR, model.SEASON, default=0)
model.ROR = Param(model.RESERVOIR, default=0)
## NEW SOS ##
model.BIGM = Param( default= 1e6 )
model.BIGM2 = Param( default= 1e5 )

###########     Param from Previous Output       ###########
if MinGlobalYear<MinYear:
    model.PrevNC = Param(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, default=0)
    model.PrevROA = Param(model.REGION, model.TIMESLICE, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, default=0)
    model.PrevTATA = Param(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, default=0)
    model.PrevSV = Param(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, default=0)
    model.PrevDSV = Param(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, default=0)
    model.PrevAP = Param(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, default=0)
    model.PrevDTEP = Param(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, default=0)
    if EmPen!=1:
        model.PrevATEM = Param(model.REGION, model.TECHNOLOGY, model.EMISSION, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, default=0.0)
        model.PrevATE = Param(model.REGION, model.TECHNOLOGY, model.EMISSION, model.FIXEDMILPYEARS, default=0.0)
        model.PrevATEPE = Param(model.REGION, model.TECHNOLOGY, model.EMISSION, model.FIXEDMILPYEARS, default=0.0)
        model.PrevATEP = Param(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, default=0.0)
        model.PrevAE = Param(model.REGION, model.EMISSION, model.FIXEDMILPYEARS,default=0.0)
    model.PrevNNT = Param(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, default=0)    
    model.PrevRRD = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevRRS = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevRLS = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevREL = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevINF = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevX = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevY = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevSUR = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevLAG = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, default=0)
    model.PrevMRE = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, default=0)
    model.PrevmRE = Param(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, default=0)

inf={}
vol_={}
su1={}
su2={}


for res in Reservoirs:
    if res[8]==0:
        inf[res[0]]=np.linspace(0,res[13],len(res[11]))
#        inf[res[0]][0]=-400
        halfscartoupper=res[12][-1]*scartoupper/100
        halfscartolower=res[12][-1]*scartolower/100
        temp2=[]
        temp3=[]
        for valore in res[12]:
            temp2.append(valore+halfscartoupper)
            temp3.append(valore-halfscartolower)
        temp2[0]=0.0
        temp3[0]=0.0
        for y in range(MinYear,MaxYear+1):
            for s in range(1,13):
                key=res[0],y,s,'SAPP'
                vol_[key]=res[11]
                vol_[key][0]=0.0
                su1[key]=temp2
                su2[key]=temp3
                
                
                 
#########			Emissions & Penalties		#############0.040
if EmPen!=1:
    model.EmissionActivityRatio = Param(model.REGION, model.TECHNOLOGY, model.EMISSION, model.MODE_OF_OPERATION, model.YEAR, default=0)
    model.EmissionsPenalty = Param( default=EmPen) #social cost of carbon 40$/tonCO2
    model.AnnualExogenousEmission = Param(model.REGION, model.EMISSION, model.YEAR, default=0)
    model.AnnualEmissionLimit = Param(model.REGION, model.EMISSION, model.YEAR, default=99999)
    model.ModelPeriodExogenousEmission = Param(model.REGION, model.EMISSION, default=0)
    model.ModelPeriodEmissionLimit = Param(model.REGION, model.EMISSION, default=99999)
          
######################
#   Model Variables  #
######################

#########		    Capacity Variables 			############# 

model.NewCapacity = Var(model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0)
model.NumberOfNewTechnologyUnits = Var(model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeIntegers, initialize=0)#ho messo 0 al posto di 0.0

#########		    Activity Variables 			#############

model.RateOfActivity = Var(model.REGION, model.TIMESLICE, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.YEAR, domain=NonNegativeReals, initialize=0.0)
model.TotalAnnualTechnologyActivityByMode = Var(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.YEAR, domain=NonNegativeReals, initialize=0.0)
#########		    Costing Variables 			#############

model.SalvageValue = Var(model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0)
model.DiscountedSalvageValue = Var(model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0)

#########			Annual Production      #############

model.AnnualProduction = Var(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.YEAR, domain=NonNegativeReals, initialize=0.0)

#########			Emissions					#############
model.DiscountedTechnologyEmissionsPenalty = Var(model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0)
if EmPen!=1:
    model.AnnualTechnologyEmissionByMode = Var(model.REGION, model.TECHNOLOGY, model.EMISSION, model.MODE_OF_OPERATION, model.YEAR, domain=NonNegativeReals, initialize=0.0)
    model.AnnualTechnologyEmission = Var(model.REGION, model.TECHNOLOGY, model.EMISSION, model.YEAR, domain=NonNegativeReals, initialize=0.0)
    model.AnnualTechnologyEmissionPenaltyByEmission = Var(model.REGION, model.TECHNOLOGY, model.EMISSION, model.YEAR, domain=NonNegativeReals, initialize=0.0)
    model.AnnualTechnologyEmissionsPenalty = Var(model.REGION, model.TECHNOLOGY, model.YEAR, domain=NonNegativeReals, initialize=0.0)
    model.AnnualEmissions = Var(model.REGION, model.EMISSION, model.YEAR, domain=NonNegativeReals, initialize=0.0)
    model.ModelPeriodEmissions = Var(model.REGION, model.EMISSION, domain=NonNegativeReals, initialize=0.0)

#########			Hydro   					  #############
model.RateOfReservoirDischarge = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.RateOfReservoirSpillage = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION,  domain=NonNegativeReals, initialize=0.0)
model.ReservoirLevelSeason = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.Release = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.Inflow = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.ReservoirLevelFinal = Var(model.RESERVOIR, domain=NonNegativeReals, initialize=0.0)
model.x = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.y = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.Surface = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.beta = Var(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.J, domain=Binary, initialize=0)
model.alpha = Var(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.I, domain=Binary, initialize=0)
model.ReleaseLag = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)
model.MaxRel = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, domain=NonNegativeReals, initialize=0.0)
model.MinRel = Var(model.RESERVOIR, model.MILPYEARS, model.SEASON, domain=NonNegativeReals, initialize=0.0)
model.ReservoirLevelFinalFineStep = Var(model.RESERVOIR, domain=NonNegativeReals, initialize=0.0)
model.CostruiteBIGM = Var(model.RESERVOIR, model.MILPYEARS, domain=Binary, initialize=0)
model.CostruiteAux = Var(model.RESERVOIR, model.MILPYEARS, domain=Binary, initialize=0)
model.gamma = Var(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, domain=NonNegativeReals, initialize=0.0)

######################
# Objective Function #
######################
model.Obj_AllYears_nolast=Var(domain=NonNegativeReals, initialize=0.0)
model.Obj_AllYearsNoXY=Var(domain=NonNegativeReals, initialize=0.0)
model.Obj_AllYearsNoXY_nolast=Var(domain=NonNegativeReals, initialize=0.0)

def ObjectiveFunction_rule(model): 
	return sum(sum(sum(((sum(model.NewCapacity[r,t,yy] for yy in model.YEAR if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0)))                      
                      + model.ResidualCapacity[r,t,y])*model.FixedCost[r,t,y]  
    + sum(sum(model.RateOfActivity[r,l,t,m,y]*model.YearSplit[l,y] for l in model.TIMESLICE)*model.VariableCost[r,t,m,y] for m in model.MODE_OF_OPERATION))
    /((1+model.DiscountRate[r])**(y-MinYearToDiscount + 0.5))
    + model.CapitalCost[r,t,y]*model.NewCapacity[r,t,y]/((1+model.DiscountRate[r])**(y-MinYearToDiscount)) 
    - model.DiscountedSalvageValue[r,t,y]+ model.DiscountedTechnologyEmissionsPenalty[r,t,y]
    for t in model.TECHNOLOGY) for y in model.YEAR) for r in model.REGION) + \
    sum(sum(sum(sum(model.x[v[0],y,ls,r] + model.y[v[0],y,ls,r] for ls in model.SEASON) for y in model.MILPYEARS) for r in model.REGION)*v[18]*0.45 for v in Reservoirs ) * model.RhoWater * model.Gravity
model.OBJ = Objective(rule=ObjectiveFunction_rule, sense=minimize)

def Obj_AllYears_nolastrule(model):
	return sum(sum(sum(((sum(model.NewCapacity[r,t,yy] for yy in model.YEARSMINUSLAST if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0)))                      
                      + model.ResidualCapacity[r,t,y])*model.FixedCost[r,t,y]  
    + sum(sum(model.RateOfActivity[r,l,t,m,y]*model.YearSplit[l,y] for l in model.TIMESLICE)*model.VariableCost[r,t,m,y] for m in model.MODE_OF_OPERATION))
    /((1+model.DiscountRate[r])**(y-MinYearToDiscount + 0.5))
    + model.CapitalCost[r,t,y]*model.NewCapacity[r,t,y]/((1+model.DiscountRate[r])**(y-MinYearToDiscount)) 
    - model.DiscountedSalvageValue[r,t,y]+ model.DiscountedTechnologyEmissionsPenalty[r,t,y]
    for t in model.TECHNOLOGY) for y in model.YEARSMINUSLAST) for r in model.REGION) + \
    sum(sum(sum(sum(model.x[v[0],y,ls,r] + model.y[v[0],y,ls,r] for ls in model.SEASON) for y in model.MILPYEARS) for r in model.REGION)*v[18]*0.45 for v in Reservoirs ) * model.RhoWater * model.Gravity == model.Obj_AllYears_nolast
model.Obj1constraint = Constraint(rule=Obj_AllYears_nolastrule) 

def Obj_AllYearsNoXYrule(model):
	return sum(sum(sum(((sum(model.NewCapacity[r,t,yy] for yy in model.YEAR if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0)))                      
                      + model.ResidualCapacity[r,t,y])*model.FixedCost[r,t,y]  
    + sum(sum(model.RateOfActivity[r,l,t,m,y]*model.YearSplit[l,y] for l in model.TIMESLICE)*model.VariableCost[r,t,m,y] for m in model.MODE_OF_OPERATION))
    /((1+model.DiscountRate[r])**(y-MinYearToDiscount + 0.5))
    + model.CapitalCost[r,t,y]*model.NewCapacity[r,t,y]/((1+model.DiscountRate[r])**(y-MinYearToDiscount)) 
    - model.DiscountedSalvageValue[r,t,y]+ model.DiscountedTechnologyEmissionsPenalty[r,t,y]
    for t in model.TECHNOLOGY) for y in model.YEAR) for r in model.REGION) == model.Obj_AllYearsNoXY
model.Obj2constraint = Constraint(rule=Obj_AllYearsNoXYrule)

def Obj_AllYearsNoXY_nolastrule(model):
	return sum(sum(sum(((sum(model.NewCapacity[r,t,yy] for yy in model.YEARSMINUSLAST if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0)))                      
                      + model.ResidualCapacity[r,t,y])*model.FixedCost[r,t,y]  
    + sum(sum(model.RateOfActivity[r,l,t,m,y]*model.YearSplit[l,y] for l in model.TIMESLICE)*model.VariableCost[r,t,m,y] for m in model.MODE_OF_OPERATION))
    /((1+model.DiscountRate[r])**(y-MinYearToDiscount + 0.5))
    + model.CapitalCost[r,t,y]*model.NewCapacity[r,t,y]/((1+model.DiscountRate[r])**(y-MinYearToDiscount)) 
    - model.DiscountedSalvageValue[r,t,y]+ model.DiscountedTechnologyEmissionsPenalty[r,t,y]
    for t in model.TECHNOLOGY) for y in model.YEARSMINUSLAST) for r in model.REGION) == model.Obj_AllYearsNoXY_nolast
model.Obj3constraint = Constraint(rule=Obj_AllYearsNoXY_nolastrule)     

#########   		Emissions Accounting		##############
if EmPen!=1:
    def AverageAnnualRateOfActivity_rule(model,r,t,m,y):
    	return sum(model.RateOfActivity[r,l,t,m,y]*model.YearSplit[l,y] for l in model.TIMESLICE) == model.TotalAnnualTechnologyActivityByMode[r,t,m,y]
    model.AverageAnnualRateOfActivity = Constraint(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.ACTIVEYEARS, rule=AverageAnnualRateOfActivity_rule)
    
    def AnnualEmissionProductionByMode_rule(model,r,t,e,m,y):
    	if model.EmissionActivityRatio[r,t,e,m,y] != 0:
    		return model.EmissionActivityRatio[r,t,e,m,y]*model.TotalAnnualTechnologyActivityByMode[r,t,m,y] == model.AnnualTechnologyEmissionByMode[r,t,e,m,y]
    	else:
    		return model.AnnualTechnologyEmissionByMode[r,t,e,m,y] == 0
    model.AnnualEmissionProductionByMode = Constraint(model.REGION, model.TECHNOLOGY, model.EMISSION, model.MODE_OF_OPERATION, model.ACTIVEYEARS, rule=AnnualEmissionProductionByMode_rule)
    
    def AnnualEmissionProduction_rule(model,r,t,e,y):
    	return sum(model.AnnualTechnologyEmissionByMode[r,t,e,m,y] for m in model.MODE_OF_OPERATION) == model.AnnualTechnologyEmission[r,t,e,y]
    model.AnnualEmissionProduction = Constraint(model.REGION, model.TECHNOLOGY, model.EMISSION, model.ACTIVEYEARS, rule=AnnualEmissionProduction_rule)
    
    def EmissionPenaltyByTechAndEmission_rule(model,r,t,e,y):
    	return model.AnnualTechnologyEmission[r,t,e,y]*model.EmissionsPenalty == model.AnnualTechnologyEmissionPenaltyByEmission[r,t,e,y]
    model.EmissionPenaltyByTechAndEmission = Constraint(model.REGION, model.TECHNOLOGY, model.EMISSION, model.ACTIVEYEARS, rule=EmissionPenaltyByTechAndEmission_rule)
    
    def EmissionsPenaltyByTechnology_rule(model,r,t,y):
    	return sum(model.AnnualTechnologyEmissionPenaltyByEmission[r,t,e,y] for e in model.EMISSION) == model.AnnualTechnologyEmissionsPenalty[r,t,y]
    model.EmissionsPenaltyByTechnology = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=EmissionsPenaltyByTechnology_rule)
    
    def DiscountedEmissionsPenaltyByTechnology_rule(model,r,t,y):
    	return model.AnnualTechnologyEmissionsPenalty[r,t,y]/((1+model.DiscountRate[r])**(y-MinYearToDiscount+0.5)) == model.DiscountedTechnologyEmissionsPenalty[r,t,y]
    model.DiscountedEmissionsPenaltyByTechnology = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=DiscountedEmissionsPenaltyByTechnology_rule)
    
    def EmissionsAccounting1_rule(model,r,e,y):
    	return sum(model.AnnualTechnologyEmission[r,t,e,y] for t in model.TECHNOLOGY) == model.AnnualEmissions[r,e,y]
    model.EmissionsAccounting1 = Constraint(model.REGION, model.EMISSION, model.ACTIVEYEARS, rule=EmissionsAccounting1_rule)
    
    def EmissionsAccounting2_rule(model,r,e):
    	return sum(model.AnnualEmissions[r,e,y] for y in model.YEAR) == model.ModelPeriodEmissions[r,e] - model.ModelPeriodExogenousEmission[r,e]
    model.EmissionsAccounting2 = Constraint(model.REGION, model.EMISSION, rule=EmissionsAccounting2_rule)
    
    def AnnualEmissionsLimit_rule(model,r,e,y):
    	return model.AnnualEmissions[r,e,y] + model.AnnualExogenousEmission[r,e,y] <= model.AnnualEmissionLimit[r,e,y]
    model.AnnualEmissionsLimit = Constraint(model.REGION, model.EMISSION, model.YEAR, rule=AnnualEmissionsLimit_rule)
    
    def ModelPeriodEmissionsLimit_rule(model,r,e):
    	return model.ModelPeriodEmissions[r,e] <= model.ModelPeriodEmissionLimit[r,e]
    model.ModelPeriodEmissionsLimit = Constraint(model.REGION, model.EMISSION, rule=ModelPeriodEmissionsLimit_rule)
else:
    def DiscountedEmissionsPenaltyByTechnology_rule(model,r,t,y):
    	return 0==model.DiscountedTechnologyEmissionsPenalty[r,t,y]
    model.DiscountedEmissionsPenaltyByTechnology = Constraint(model.REGION, model.TECHNOLOGY, model.YEAR, rule=DiscountedEmissionsPenaltyByTechnology_rule)    

#####################
# Constraints       #
#####################

##########      Annual Production            ##############

def AnnualProduction_rule(model,r,t,m,y):
    return sum(sum(model.RateOfActivity[r,l,t,m,y]*model.OutputActivityRatio[r,t,f,m,y] for f in model.FUEL)*model.YearSplit[l,y] for l in model.TIMESLICE ) == model.AnnualProduction[r,t,m,y]
model.AnnualProduction_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.ACTIVEYEARS, rule=AnnualProduction_rule)

#########       	Capacity Adequacy A	     	#############
def ConstraintCapacity_rule(model,r,l,t,y):
    return sum(model.RateOfActivity[r,l,t,m,y] for m in model.MODE_OF_OPERATION) <= (sum(model.NewCapacity[r,t,yy] for yy in model.YEAR if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0))) + model.ResidualCapacity[r,t,y])*model.CapacityFactor[r,t,l,y]*model.CapacityToActivityUnit[r,t]
model.ConstraintCapacity = Constraint(model.REGION, model.TIMESLICE, model.TECHNOLOGY, model.ACTIVEYEARS, rule=ConstraintCapacity_rule)

#########       	Capacity Adequacy B		 	#############
def PlannedMaintenance_rule(model,r,t,y):
	 return sum(sum(model.RateOfActivity[r,l,t,m,y] for m in model.MODE_OF_OPERATION)*model.YearSplit[l,y] for l in model.TIMESLICE) <= sum((sum(model.NewCapacity[r,t,yy] for yy in model.YEAR if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0))) + model.ResidualCapacity[r,t,y])*model.CapacityFactor[r,t,l,y]*model.YearSplit[l,y] for l in model.TIMESLICE)*model.AvailabilityFactor[r,t,y]*model.CapacityToActivityUnit[r,t]

model.PlannedMaintenance = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=PlannedMaintenance_rule)

#########	        Energy Balance A    	 	#############

def EnergyBalanceEachTS5_rule(model,r,l,f,y):
    return sum(sum(model.RateOfActivity[r,l,t,m,y]*model.OutputActivityRatio[r,t,f,m,y] for m in model.MODE_OF_OPERATION ) for t in model.TECHNOLOGY)*model.YearSplit[l,y] >= model.SpecifiedAnnualDemand[r,f,y]*model.SpecifiedDemandProfile[r,f,l,y] + sum(sum(model.RateOfActivity[r,l,t,m,y]*model.InputActivityRatio[r,t,f,m,y] for m in model.MODE_OF_OPERATION ) for t in model.TECHNOLOGY)*model.YearSplit[l,y]
model.EnergyBalanceEachTS5 = Constraint(model.REGION, model.TIMESLICE, model.FUEL, model.ACTIVEYEARS, rule=EnergyBalanceEachTS5_rule)

#########        	Energy Balance B		 	#############

def EnergyBalanceEachYear4_rule(model,r,f,y):
	 return sum(sum(sum(model.RateOfActivity[r,l,t,m,y]*model.OutputActivityRatio[r,t,f,m,y] for m in model.MODE_OF_OPERATION ) for t in model.TECHNOLOGY)*model.YearSplit[l,y] for l in model.TIMESLICE) >= sum(sum(sum(model.RateOfActivity[r,l,t,m,y]*model.InputActivityRatio[r,t,f,m,y] for m in model.MODE_OF_OPERATION ) for t in model.TECHNOLOGY)*model.YearSplit[l,y] for l in model.TIMESLICE)
model.EnergyBalanceEachYear4 = Constraint(model.REGION, model.FUEL, model.ACTIVEYEARS, rule=EnergyBalanceEachYear4_rule)

#########      		Total Capacity Constraints 	##############
def TotalAnnualMaxCapacityConstraint_rule(model,r,t,y): 
	return sum(model.NewCapacity[r,t,yy] for yy in model.YEAR if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0))) + model.ResidualCapacity[r,t,y] <= model.TotalAnnualMaxCapacity[r,t,y]
model.TotalAnnualMaxCapacityConstraint = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=TotalAnnualMaxCapacityConstraint_rule)

def TotalAnnualMinCapacityConstraint_rule(model,r,t,y):
    return sum(model.NewCapacity[r,t,yy] for yy in model.YEAR if ((y-yy < model.OperationalLife[r,t]) and (y-yy >= 0))) + model.ResidualCapacity[r,t,y] >= model.TotalAnnualMinCapacity[r,t,y]

model.TotalAnnualMinCapacityConstraint = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=TotalAnnualMinCapacityConstraint_rule)           

def TotalNewCapacity_2_rule(model,r,t,y):
    if model.CapacityOfOneTechnologyUnit[r,t,y] != 0.0: #parametro definito solo per le nuove centrali, pari alla capacity prevista per la nuova centrale
        return model.CapacityOfOneTechnologyUnit[r,t,y]*model.NumberOfNewTechnologyUnits[r,t,y] == model.NewCapacity[r,t,y] #la variabile dunque sarà l'anno di costruzione
    else: 
        return Constraint.Skip
model.TotalNewCapacity_2 = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=TotalNewCapacity_2_rule)

if VecchieOrNuove=="Nuove":
    def NNTDG(model,r,y):
        return model.NumberOfNewTechnologyUnits[r,'DGU',y]==model.NumberOfNewTechnologyUnits[r,'DGL',y]
    model.NNTDG_Rule = Constraint(model.REGION, model.ACTIVEYEARS, rule=NNTDG)
    
    def NNTBG(model,r,y):
        return model.NumberOfNewTechnologyUnits[r,'BGN',y]==model.NumberOfNewTechnologyUnits[r,'BGS',y]
    model.NNTBG_Rule = Constraint(model.REGION, model.ACTIVEYEARS, rule=NNTBG)
    
    def ANNUALDG(model,r,m,y):
        return model.AnnualProduction[r,'DGU',m,y]==model.AnnualProduction[r,'DGL',m,y]
    model.ANNUALDG_Rule = Constraint(model.REGION, model.MODE_OF_OPERATION, model.ACTIVEYEARS, rule=ANNUALDG)
    
    def ANNUALBG(model,r,m,y):
        return model.AnnualProduction[r,'BGS',m,y]==model.AnnualProduction[r,'BGN',m,y]
    model.ANNUALBG_Rule = Constraint(model.REGION, model.MODE_OF_OPERATION, model.ACTIVEYEARS, rule=ANNUALBG)


#########           Salvage Value            	#############

def SalvageValueAtEndOfPeriod1_rule(model,r,t,y): 
	if model.DepreciationMethod[r] == 1 and ((y + model.OperationalLife[r,t]-1) > MaxYearToDiscount) and model.DiscountRate[r]>0: 
		return model.SalvageValue[r,t,y] == model.CapitalCost[r,t,y]*model.NewCapacity[r,t,y]*(1-(((1+model.DiscountRate[r])**(MaxYearToDiscount- y+1)-1)/((1+model.DiscountRate[r])**model.OperationalLife[r,t]-1)))
	elif (model.DepreciationMethod[r] == 1 and ((y + model.OperationalLife[r,t]-1) > MaxYearToDiscount) and model.DiscountRate[r] == 0) or (model.DepreciationMethod[r] == 2 and (y + model.OperationalLife[r,t]-1) > (MaxYearToDiscount)):
		return model.SalvageValue[r,t,y] == model.CapitalCost[r,t,y]*model.NewCapacity[r,t,y]*(1-(MaxYearToDiscount- y+1)/model.OperationalLife[r,t])
	else:
		return model.SalvageValue[r,t,y] == 0
   
model.SalvageValueAtEndOfPeriod1 = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=SalvageValueAtEndOfPeriod1_rule)

def SalvageValueDiscountedToStartYear_rule(model,r,t,y):
    return model.DiscountedSalvageValue[r,t,y] == model.SalvageValue[r,t,y]/((1+model.DiscountRate[r])**(1+MaxYearToDiscount-MinYearToDiscount))
model.SalvageValueDiscountedToStartYear = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=SalvageValueDiscountedToStartYear_rule)

#########    		New Capacity Constraints  	##############

def TotalAnnualMaxNewCapacityConstraint_rule(model,r,t,y):
    return model.NewCapacity[r,t,y] <= model.TotalAnnualMaxCapacityInvestment[r,t,y]
model.TotalAnnualMaxNewCapacityConstraint = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=TotalAnnualMaxNewCapacityConstraint_rule)

def TotalAnnualMinNewCapacityConstraint_rule(model,r,t,y):
    return model.NewCapacity[r,t,y] >= model.TotalAnnualMinCapacityInvestment[r,t,y]
model.TotalAnnualMinNewCapacityConstraint = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=TotalAnnualMinNewCapacityConstraint_rule)       

#########   		Annual Activity Constraints	##############

def TotalAnnualTechnologyActivityUpperLimit_rule(model,r,t,y):
    return sum(sum(model.RateOfActivity[r,l,t,m,y] for m in model.MODE_OF_OPERATION)*model.YearSplit[l,y] for l in model.TIMESLICE) <= model.TotalTechnologyAnnualActivityUpperLimit[r,t,y]
model.TotalAnnualTechnologyActivityUpperlimit = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=TotalAnnualTechnologyActivityUpperLimit_rule)

def TotalAnnualTechnologyActivityLowerLimit_rule(model,r,t,y):
#    if model.TotalTechnologyAnnualActivityLowerLimit[r,t,y] != 0:
    return sum(sum(model.RateOfActivity[r,l,t,m,y] for m in model.MODE_OF_OPERATION)*model.YearSplit[l,y] for l in model.TIMESLICE) >= model.TotalTechnologyAnnualActivityLowerLimit[r,t,y]
#    else:
#        return Constraint.Skip
model.TotalAnnualTechnologyActivityLowerlimit = Constraint(model.REGION, model.TECHNOLOGY, model.ACTIVEYEARS, rule=TotalAnnualTechnologyActivityLowerLimit_rule)

#########    		Total Activity Constraints 	##############

def TotalModelHorizonTechnologyActivityUpperLimit_rule(model,r,t):
    if model.TotalTechnologyModelPeriodActivityUpperLimit[r,t] != 0:
        return sum(sum(sum(model.RateOfActivity[r,l,t,m,y] for m in model.MODE_OF_OPERATION)*model.YearSplit[l,y] for l in model.TIMESLICE) for y in model.YEAR) <= model.TotalTechnologyModelPeriodActivityUpperLimit[r,t]
    else:
        return Constraint.Skip
model.TotalModelHorizonTechnologyActivityUpperLimit = Constraint(model.REGION, model.TECHNOLOGY, rule=TotalModelHorizonTechnologyActivityUpperLimit_rule)

def TotalModelHorizonTechnologyActivityLowerLimit_rule(model,r,t):
    if model.TotalTechnologyModelPeriodActivityLowerLimit[r,t] != 0:
        return sum(sum(sum(model.RateOfActivity[r,l,t,m,y] for m in model.MODE_OF_OPERATION)*model.YearSplit[l,y] for l in model.TIMESLICE) for y in model.YEAR) >= model.TotalTechnologyModelPeriodActivityLowerLimit[r,t]
    else:
        return Constraint.Skip
model.TotalModelHorizonTechnologyActivityLowerLimit = Constraint(model.REGION, model.TECHNOLOGY, rule=TotalModelHorizonTechnologyActivityLowerLimit_rule)

#########   		Hydro Equations	##############

def Hy1_RateOfReservoirDischarge_rule(model,v,y,ls,r):
    return (sum(sum(sum(model.RateOfActivity[r,l,t,m,y] * model.TechnologyFromHydro[r,t,m,v] / model.HydroTurbineEfficiency[t] for m 
                       in model.MODE_OF_OPERATION if (model.TechnologyFromHydro[r,t,m,v] != 0))
    *model.Conversionls[l,ls] for l in model.TIMESLICE)
    /model.RhoWater/model.Gravity
    /model.WattsToModelUnits
    /model.CapacityToActivityUnit[r,t] for t in model.TECHNOLOGY)
    /model.ReservoirHead[v]) == model.RateOfReservoirDischarge[v,y,ls,r]
model.Hy1_RateOfReservoirDischarge = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Hy1_RateOfReservoirDischarge_rule)

def Hy9_and_Hy10_ReservoirLevelSeason_rule(model,v,y,ls,r): #controlla cosa sono quegli 1e-3 della riga 367-->è il /1000 che considera evap in mm/m2
    if sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0 and sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0:
        return Constraint.Skip
    else:
        if (ls==min(model.SEASON) and y==min(model.YEAR)):
            return model.ReservoirLevelStart[r,v] == model.ReservoirLevelSeason[v,y,ls,r]
        elif (ls==min(model.SEASON) and y!=min(model.YEAR)):
            return model.ReservoirLevelSeason[v,y-1,max(model.SEASON),r] - 1e-3 * model.Surface[v,y-1,max(model.SEASON),r] * model.Evap[v,max(model.SEASON)] + 1e-9 * sum((model.Inflow[v,y-1,max(model.SEASON),r] - model.Release[v,y-1,max(model.SEASON),r]) * model.FlowUnittoYearConversion * model.YearSplit[l,y-1] * model.Conversionls[l,max(model.SEASON)] for l in model.TIMESLICE if model.Conversionls[l,max(model.SEASON)]>0) == model.ReservoirLevelSeason[v,y,ls,r]
        else:
            return model.ReservoirLevelSeason[v,y,ls-1,r] - 1e-3 * model.Surface[v,y,ls-1,r] * model.Evap[v,ls-1] + 1e-9 * sum((model.Inflow[v,y,ls-1,r] - model.Release[v,y,ls-1,r]) * model.FlowUnittoYearConversion * model.YearSplit[l,y] * model.Conversionls[l,ls-1] for l in model.TIMESLICE if model.Conversionls[l,ls-1]>0) == model.ReservoirLevelSeason[v,y,ls,r]
model.Hy9_and_Hy10_ReservoirLevelSeason = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Hy9_and_Hy10_ReservoirLevelSeason_rule)


def Kariba_rule1(model,y,ls):
   return sum(sum(sum(model.RateOfActivity[r,l,'KAS',m,y] * model.TechnologyFromHydro[r,'KAS',m,'RESKA'] / model.HydroTurbineEfficiency['KAS'] for m 
                      in model.MODE_OF_OPERATION if (model.TechnologyFromHydro[r,'KAS',m,'RESKA'] != 0)) * model.Conversionls[l,ls] for l in model.TIMESLICE) 
   /model.RhoWater/model.Gravity 
   /model.WattsToModelUnits 
   /model.CapacityToActivityUnit[r,'KAS'] for r in model.REGION) \
   /model.ReservoirHead['RESKA'] >= model.RateOfReservoirDischarge['RESKA',y,ls,'SAPP'] * 0.512
model.KaribaPolicyConstraint1 = Constraint(model.ACTIVEMILPYEARS, model.SEASON, rule=Kariba_rule1)

def Kariba_rule2(model,y,ls):
   return sum(sum(sum(model.RateOfActivity[r,l,'KAS',m,y] * model.TechnologyFromHydro[r,'KAS',m,'RESKA'] / model.HydroTurbineEfficiency['KAS'] for m 
                      in model.MODE_OF_OPERATION if (model.TechnologyFromHydro[r,'KAS',m,'RESKA'] != 0)) * model.Conversionls[l,ls] for l in model.TIMESLICE) 
   /model.RhoWater/model.Gravity 
   /model.WattsToModelUnits 
   /model.CapacityToActivityUnit[r,'KAS'] for r in model.REGION) \
   /model.ReservoirHead['RESKA'] <= model.RateOfReservoirDischarge['RESKA',y,ls,'SAPP'] * 0.556
model.KaribaPolicyConstraint2 = Constraint(model.ACTIVEMILPYEARS, model.SEASON, rule=Kariba_rule2)

#########   		Hydro Constraint	##############

def MaxStorage_rule(model,v,y,ls,r):
    return model.ReservoirLevelSeason[v,y,ls,r] <= model.ReservoirStorageVolume[v]
model.MaxStorage = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=MaxStorage_rule)

def FinalState_rule(model,v,r):
    return model.ReservoirLevelSeason[v,max(model.ACTIVEMILPYEARS),max(model.SEASON),r] - 1e-3 * model.Surface[v,max(model.ACTIVEMILPYEARS),max(model.SEASON),r] * model.Evap[v,max(model.SEASON)] + 1e-9 * sum((model.Inflow[v,max(model.ACTIVEMILPYEARS),max(model.SEASON),r] - model.Release[v,max(model.ACTIVEMILPYEARS),max(model.SEASON),r]) * model.FlowUnittoYearConversion * model.YearSplit[l,max(model.ACTIVEMILPYEARS)] * model.Conversionls[l,max(model.SEASON)] for l in model.TIMESLICE if model.Conversionls[l,max(model.SEASON)]>0) == model.ReservoirLevelFinal[v]
model.FinalState = Constraint(model.RESERVOIR, model.REGION, rule=FinalState_rule)

def FinalState2_rule(model,v):
    return model.ReservoirLevelFinal[v] <= model.ReservoirStorageVolume[v]
model.FinalState2 = Constraint(model.RESERVOIR, rule=FinalState2_rule)

def FineStepState_rule(model,v,r):
    return model.ReservoirLevelSeason[v,(max(model.ACTIVEMILPYEARS)-OverlappingSize),max(model.SEASON),r] - 1e-3 * model.Surface[v,(max(model.ACTIVEMILPYEARS)-OverlappingSize),max(model.SEASON),r] * model.Evap[v,max(model.SEASON)] + 1e-9 * sum((model.Inflow[v,(int(max(model.ACTIVEMILPYEARS))-OverlappingSize),max(model.SEASON),r] - model.Release[v,(int(max(model.ACTIVEMILPYEARS))-OverlappingSize),max(model.SEASON),r]) * model.FlowUnittoYearConversion * model.YearSplit[l,(int(max(model.ACTIVEMILPYEARS))-OverlappingSize)] * model.Conversionls[l,max(model.SEASON)] for l in model.TIMESLICE if model.Conversionls[l,max(model.SEASON)]>0) == model.ReservoirLevelFinalFineStep[v]
model.FineStepState = Constraint(model.RESERVOIR, model.REGION, rule=FineStepState_rule)

def FineStepState2_rule(model,v):
    return model.ReservoirLevelFinalFineStep[v] <= model.ReservoirStorageVolume[v]
model.FineStepState2 = Constraint(model.RESERVOIR, rule=FineStepState2_rule)

def nCostruite(model,v,y,r):
    return sum(sum(model.NumberOfNewTechnologyUnits[r,t,yy] for yy in model.MILPYEARS if (y-yy >= 0)) for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0)) == model.CostruiteAux[v,y] + model.CostruiteBIGM[v,y]
model.nCostruiteConstraint = Constraint(model.RESERVOIR, model.MILPYEARS, model.REGION, rule=nCostruite)    

def CostruiteXY(model,v,y):
    return model.CostruiteBIGM[v,y] >= model.CostruiteAux[v,y]
model.CostruiteXYConstraint = Constraint(model.RESERVOIR, model.MILPYEARS, rule=CostruiteXY)    

#########   		Hydro Discharge requirements ##############

def ReleaseLag_rule(model,vv,y,ls,r,v): #vv=monte v=valle 
    if model.DownstreamReservoirTag[v,vv]!=0:
        if ls in range(1,value(model.InflowLagTime[v,vv])+1):
            if y==min(model.ACTIVEMILPYEARS):
                return model.StartRelease[vv,max(model.SEASON)-model.InflowLagTime[v,vv]+ls] == model.ReleaseLag[vv,y,ls,r] #appena aggiunto il v finale
            else:
                return model.Release[vv,y-1,max(model.SEASON)-model.InflowLagTime[v,vv]+ls,r] == model.ReleaseLag[vv,y,ls,r] #appena aggiunto il v finale
        else:
             return model.Release[vv,y,ls-model.InflowLagTime[v,vv],r] == model.ReleaseLag[vv,y,ls,r]
    else:
        return Constraint.Skip
model.ReleaseLagConstraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, model.RESERVOIR, rule=ReleaseLag_rule)
     
def Inflow_rule(model,v,y,ls,r):
    return model.ReservoirExternalInflow[r,ls,y,v] + sum(model.ReleaseLag[vv,y,ls,r]*model.DownstreamReservoirTag[v,vv] for vv in model.RESERVOIR) == model.Inflow[v,y,ls,r]
model.InflowConstraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Inflow_rule)

def Release_rule(model,v,y,ls,r):
    return model.RateOfReservoirDischarge[v,y,ls,r] + model.RateOfReservoirSpillage[v,y,ls,r] == model.Release[v,y,ls,r]                
model.ReleaseConstraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Release_rule)

def Spillage_rule(model,v,y,ls,r): #per le centrali nuove nel periodo dove non posso essere costruite
    if sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0 and sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0:
        return model.Inflow[v,y,ls,r] == model.RateOfReservoirSpillage[v,y,ls,r]
    else:
        return Constraint.Skip
model.Spillage_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Spillage_rule)

def Release_rule_1(model,v,y,ls,r): #primo vincolo release centrali costruibili in questo run
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.RateOfReservoirSpillage[v,y,ls,r] <= model.Inflow[v,y,ls,r] + model.BIGM2 * model.CostruiteBIGM[v,y]
    else:
        return Constraint.Skip
model.ReleaseConstraint1 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Release_rule_1)

def Release_rule_4(model,v,y,ls,r): #primo vincolo release centrali costruibili in questo run
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.RateOfReservoirSpillage[v,y,ls,r] >= model.Inflow[v,y,ls,r] - model.BIGM2 * model.CostruiteBIGM[v,y]
    else:
        return Constraint.Skip
model.ReleaseConstraint4 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Release_rule_4)

def Release_rule_2(model,v,y,ls,r): #primo vincolo release centrali costruibili in questo run
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.RateOfReservoirSpillage[v,y,ls,r] <= model.Release[v,y,ls,r] - model.RateOfReservoirDischarge[v,y,ls,r] + model.BIGM2 * (1-model.CostruiteBIGM[v,y])
    else:
        return Constraint.Skip
model.ReleaseConstraint2 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Release_rule_2)
   
def Release_rule_3(model,v,y,ls,r): #primo vincolo release centrali costruibili in questo run
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.RateOfReservoirSpillage[v,y,ls,r] >= model.Release[v,y,ls,r] - model.RateOfReservoirDischarge[v,y,ls,r] - model.BIGM2 * (1-model.CostruiteBIGM[v,y])
    else:
        return Constraint.Skip
model.ReleaseConstraint3 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=Release_rule_3)        
    
def MEF_rule(model,v,y,ls,r): #N.B: ho valori di MEF solo per VF e ITT
    if model.ROR[v]==1: #se è un ROR NON HO STORAGE! quello che entra=auello che esce
        if model.MEF[v,ls] < value(model.ReservoirExternalInflow[r,ls,y,v]): #se quello che entra è > del MEF che deve uscire 
            return model.RateOfReservoirSpillage[v,y,ls,r] >= model.MEF[v,ls] #impongo che esca una quantità ALMENO PARI al MEF
        else: #se il MEF è maggiore di quello che entra
            return model.RateOfReservoirSpillage[v,y,ls,r] >= model.ReservoirExternalInflow[r,ls,y,v] 
            #impongo che deve uscire una quantità ALMENO PARI A QUELLA CHE ENTRA
    elif model.ROR[v]==0: #se la nostra diga HA STORAGE
        return model.Release[v,y,ls,r] >= model.MEF[v,ls]
    else:
        return Constraint.Skip
model.MEFConstraint = Constraint(model.RESERVOIR,model.ACTIVEMILPYEARS,model.SEASON,model.REGION,rule=MEF_rule)

def x_nuoveperiodo1019(model,v,y,ls,r):
    if sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0 and sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0:
        return 0 == model.x[v,y,ls,r]
    else:
        return Constraint.Skip
model.x_periodo = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=x_nuoveperiodo1019)
 
def y_nuoveperiodo1019(model,v,y,ls,r):
    if sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0 and sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0:
        return 0 == model.y[v,y,ls,r]
    else:
        return Constraint.Skip
model.y_periodo = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=y_nuoveperiodo1019)
 
def xyC_rule(model,v,y,ls,r): #vincolo rule curve centrali vecchie o nuove se già costruite in un run precedente
    if sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0.0: #caso centrali costruite, sia preesistenti sia nuoves
        return model.ReservoirLevelSeason[v,y,ls,r] - model.RuleCurve[v,y,ls] == model.x[v,y,ls,r] - model.y[v,y,ls,r]
    elif sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0: #centrali nuove dal 2010, al 2018, quando non possono essere costruite
        return model.ReservoirLevelSeason[v,y,ls,r] == model.ReservoirLevelStart[r,v] #in teoria così impongo x e y a 0 nel periodo dove non posso costruire le centrali nuove
    else:
        return Constraint.Skip
model.xyC = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xyC_rule)
       
#vincoli rule curve centrali nuove dall'anno di costruzione in poi
    
def xC_rule_nuove1(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.x[v,y,ls,r]<= model.BIGM2*model.CostruiteBIGM[v,y]
    else:
        return Constraint.Skip
model.xC_nuove1 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xC_rule_nuove1)

def xC_rule_nuove2(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.x[v,y,ls,r]>= -model.BIGM2*model.CostruiteBIGM[v,y]
    else:
        return Constraint.Skip
model.xC_nuove2 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xC_rule_nuove2)         

def yC_rule_nuove1(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.y[v,y,ls,r]<= model.BIGM2*model.CostruiteBIGM[v,y]
    else:
        return Constraint.Skip
model.yC_nuove1 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=yC_rule_nuove1)

def yC_rule_nuove2(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.y[v,y,ls,r]>= -model.BIGM2*model.CostruiteBIGM[v,y]
    else:
        return Constraint.Skip
model.yC_nuove2 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=yC_rule_nuove2)         

def xyC_rule_nuove1(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.ReservoirLevelSeason[v,y,ls,r] - model.ReservoirLevelStart[r,v] <= model.BIGM2*(model.CostruiteBIGM[v,y])
    else:
        return Constraint.Skip
model.xyC_nuove1 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xyC_rule_nuove1)

def xyC_rule_nuove2(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.ReservoirLevelSeason[v,y,ls,r] - model.ReservoirLevelStart[r,v] >= -model.BIGM2*(model.CostruiteBIGM[v,y])
    else:
        return Constraint.Skip
model.xyC_nuove2 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xyC_rule_nuove2)     

def xyC_rule_nuove3(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.x[v,y,ls,r] - model.y[v,y,ls,r] <=(model.ReservoirLevelSeason[v,y,ls,r] - model.RuleCurve[v,y,ls]) + model.BIGM2*(1-model.CostruiteBIGM[v,y])
    else:
        return Constraint.Skip
model.xyC_nuove3 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xyC_rule_nuove3)

def xyC_rule_nuove4(model,v,y,ls,r): #vincolo rule curve centrali nuove
    if (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):  #centrali preesistenti e centrali nuove SE NON GIA' COSTRUITE nei run precedenti
        return model.x[v,y,ls,r] - model.y[v,y,ls,r] >= (model.ReservoirLevelSeason[v,y,ls,r] - model.RuleCurve[v,y,ls]) - model.BIGM2*(1-model.CostruiteBIGM[v,y])
    else:
        return Constraint.Skip
model.xyC_nuove4 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xyC_rule_nuove4)     


def xyC25_rule(model,v,y,ls,r): #vincolo rule curve centrali vecchie o nuove se già costruite in un run precedente
    return model.ReservoirLevelSeason[v,y,ls,r] >= model.RuleCurve[v,y,ls]/2
model.xyC25 = Constraint(model.OLDRESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=xyC25_rule)
 
def MinimumOperatingLevel1_rule(model,v,y,ls,r): #vincolo rule curve centrali vecchie o nuove se già costruite in un run precedente
    return model.ReservoirLevelSeason[v,y,ls,r] >= model.ReservoirStorageVolume[v]*0.1 - model.BIGM2*(1-model.gamma[v,y,ls,r])
model.MinimumOperatingLevel1 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=MinimumOperatingLevel1_rule)

def MinimumOperatingLevel2_rule(model,v,y,ls,r): #vincolo rule curve centrali vecchie o nuove se già costruite in un run precedente
    return model.ReservoirLevelSeason[v,y,ls,r] <= model.ReservoirStorageVolume[v]*0.1 + model.BIGM2*model.gamma[v,y,ls,r]
model.MinimumOperatingLevel2 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=MinimumOperatingLevel2_rule)

def MinimumOperatingLevel3_rule(model,v,y,ls,r): #vincolo rule curve centrali vecchie o nuove se già costruite in un run precedente
    return model.RateOfReservoirDischarge[v,y,ls,r] >= - model.BIGM2*model.gamma[v,y,ls,r]
model.MinimumOperatingLevel3 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=MinimumOperatingLevel3_rule)

def MinimumOperatingLevel4_rule(model,v,y,ls,r): #vincolo rule curve centrali vecchie o nuove se già costruite in un run precedente
    return model.RateOfReservoirDischarge[v,y,ls,r] <= model.BIGM2*model.gamma[v,y,ls,r]
model.MinimumOperatingLevel4 = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=MinimumOperatingLevel4_rule)


#########   Max/Min Release function SOS data and constraints ###########
vettoreNomiReservoir=[]
for res in Reservoirs:
    if res[8]==0:
        vettoreNomiReservoir.append(res[0])
        
model.HyC7_SurfaceConstraintLower = Piecewise(vettoreNomiReservoir, model.ACTIVEMILPYEARS, model.SEASON, model.REGION,
															model.Surface,
															model.ReservoirLevelSeason,
															pw_pts=vol_,
															pw_constr_type='LB',
															f_rule=su2,
															unbounded_domain_var=True,
															warn_domain_coverage=False)

model.HyC7_SurfaceConstraintUpper = Piecewise(vettoreNomiReservoir, model.ACTIVEMILPYEARS, model.SEASON, model.REGION,
															model.Surface,
															model.ReservoirLevelSeason,
															pw_pts=vol_,
															pw_constr_type='UB',
															f_rule=su1,
															unbounded_domain_var=True,
															warn_domain_coverage=False)


def beta_j_1_rule(model,v,y,ls,r):
    if (model.ROR[v]==0):
        return model.Inflow[v,y,ls,r]<= sum(model.beta[v,y,ls,j] * inf[v][int(j)] for j in model.J)
    else:
        return Constraint.Skip
model.beta_j_1_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=beta_j_1_rule)

def beta_j_2_rule(model,v,y,ls,r):
    if (model.ROR[v]==0):
        return model.Inflow[v,y,ls,r]>= sum(model.beta[v,y,ls,j] * inf[v][int(j)-1] for j in model.J)
    else:
        return Constraint.Skip
model.beta_j_2_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=beta_j_2_rule)


def beta_j_3_rule(model,v,y,ls):
    if (model.ROR[v]==0):
        return 1 == sum(model.beta[v,y,ls,j] for j in model.J)
    else:
        return Constraint.Skip
model.beta_j_3_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, rule=beta_j_3_rule)

vol={}
for res in Reservoirs:
    if res[8]==0:
        vol[res[0]]=res[11]
        
def alpha_i_1_rule(model,v,y,ls,r):
    if (model.ROR[v]==0):
        return model.ReservoirLevelSeason[v,y,ls,r]<= sum(model.alpha[v,y,ls,i] * vol[v][int(i)] for i in model.I)
    else:
        return Constraint.Skip
model.alpha_i_1_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=alpha_i_1_rule)

def alpha_i_2_rule(model,v,y,ls,r):
    if (model.ROR[v]==0):
        return model.ReservoirLevelSeason[v,y,ls,r]>= sum(model.alpha[v,y,ls,i] * vol[v][int(i)-1] for i in model.I)
    else:
        return Constraint.Skip
model.alpha_i_2_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=alpha_i_2_rule)


def alpha_i_3_rule(model,v,y,ls):
    if (model.ROR[v]==0):
        return 1 == sum(model.alpha[v,y,ls,i] for i in model.I)
    else:
        return Constraint.Skip
model.alpha_i_3_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, rule=alpha_i_3_rule)

max_rel={}

for res in Reservoirs:
    if res[8]==0:
        with open (path + res[14]) as f: #da cambiare per le nuove reservoir!!!!
            reader = csv.reader(f,delimiter=',')
            temp = enumerate(list(reader))
        temp=([el[1] for el in temp])
        max_rel[res[0]]=temp


def MaxRel_1_rule(model,v,y,ls,j):
    if (model.ROR[v]==0):
        return model.MaxRel[v,y,ls] <= sum(model.alpha[v,y,ls,i] * float(max_rel[v][i][j]) for i in model.I) + model.BIGM * (1-model.beta[v,y,ls,j])
    else:
        return Constraint.Skip
model.MaxRel_1_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.J, rule=MaxRel_1_rule)

def MaxRel_2_rule(model,v,y,ls,j):
    if (model.ROR[v]==0):
        return model.MaxRel[v,y,ls] >= sum(model.alpha[v,y,ls,i] * float(max_rel[v][i][j]) for i in model.I)- model.BIGM * (1-model.beta[v,y,ls,j])
    else:
        return Constraint.Skip
model.MaxRel_2_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.J, rule=MaxRel_2_rule)


min_rel={}

for res in Reservoirs:
    if res[8]==0:
        with open (path + res[15]) as f: #da cambiare per le nuove reservoir!!!!
            reader = csv.reader(f,delimiter=',')
            temp = enumerate(list(reader))
        temp=([el[1] for el in temp])
        min_rel[res[0]]=temp


def MinRel_1_rule(model,v,y,ls,j):
    if (model.ROR[v]==0):
        return model.MinRel[v,y,ls] <= sum(model.alpha[v,y,ls,i] * float(min_rel[v][i-1][j-1]) for i in model.I)+ model.BIGM * (1-model.beta[v,y,ls,j])
    else:
        return Constraint.Skip
model.MinRel_1_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.J, rule=MinRel_1_rule)

def MinRel_2_rule(model,v,y,ls,j):
    if (model.ROR[v]==0):
        return model.MinRel[v,y,ls] >= sum(model.alpha[v,y,ls,i] * float(min_rel[v][i-1][j-1]) for i in model.I)- model.BIGM * (1-model.beta[v,y,ls,j])
    else:
        return Constraint.Skip
model.MinRel_2_Constraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.J, rule=MinRel_2_rule)

def MaxRel_rule(model,v,y,ls,r):
    if (model.ROR[v]==0):
        #già costruite nei run prima o già presenti
        if sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0.0:
            return model.Release[v,y,ls,r] <= model.MaxRel[v,y,ls]
        #puo' costruirle nel presente run
        elif (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):
            return model.Release[v,y,ls,r] <= model.MaxRel[v,y,ls] + model.BIGM2*(1-model.CostruiteBIGM[v,y])
        else: #non puo' costruirle nel presente run
            return Constraint.Skip
    else:
        return Constraint.Skip
model.MaxRelConstraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=MaxRel_rule)

def MinRel_rule(model,v,y,ls,r):
    if (model.ROR[v]==0):
        if sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0.0:
            return model.Release[v,y,ls,r] >= model.MinRel[v,y,ls]
        elif (sum(model.CapacityOfOneTechnologyUnit[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))>0 and sum(model.ResidualCapacity[r,t,y] for t in model.TECHNOLOGY if (model.TechnologyFromHydro[r,t,1,v] != 0))==0.0):
            return model.Release[v,y,ls,r] >= model.MinRel[v,y,ls] - model.BIGM2*(1-model.CostruiteBIGM[v,y])
        else:
            return Constraint.Skip
    else:
        return Constraint.Skip
model.MinRelConstraint = Constraint(model.RESERVOIR, model.ACTIVEMILPYEARS, model.SEASON, model.REGION, rule=MinRel_rule)
    
   

if MinGlobalYear<MinYear:
    def FixedAnnualProduction_rule(model,r,t,m,y):
        return model.PrevAP[r,t,m,y] == model.AnnualProduction[r,t,m,y]
    model.FixedAnnualProduction_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, rule=FixedAnnualProduction_rule)

    def FixedTATA_rule(model,r,t,m,y):
        return model.PrevTATA[r,t,m,y] == model.TotalAnnualTechnologyActivityByMode[r,t,m,y]
    model.FixedTATA_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, rule=FixedTATA_rule)

    if EmPen!=1:
        def FixedDTEP_rule(model,r,t,y):
            return model.PrevDTEP[r,t,y] == model.DiscountedTechnologyEmissionsPenalty[r,t,y]
        model.FixedDTEP_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, rule=FixedDTEP_rule)

        def FixedATEM_rule(model,r,t,e,m,y):
            return model.PrevATEM[r,t,e,m,y] == model.AnnualTechnologyEmissionByMode[r,t,e,m,y]
        model.FixedATEM_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.EMISSION, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, rule=FixedATEM_rule)

        def FixedATE_rule(model,r,t,e,y):
            return model.PrevATE[r,t,e,y] == model.AnnualTechnologyEmission[r,t,e,y]
        model.FixedATE_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.EMISSION, model.FIXEDMILPYEARS, rule=FixedATE_rule)

        def FixedATEPE_rule(model,r,t,e,y):
            return model.PrevATEPE[r,t,e,y] == model.AnnualTechnologyEmissionPenaltyByEmission[r,t,e,y]
        model.FixedATEPE_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.EMISSION, model.FIXEDMILPYEARS, rule=FixedATEPE_rule)
        
        def FixedATEP_rule(model,r,t,y):
            return model.PrevATEP[r,t,y] == model.AnnualTechnologyEmissionsPenalty[r,t,y]
        model.FixedATEP_constraint = Constraint(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, rule=FixedATEP_rule)

        def FixedAE_rule(model,r,e,y):
            return model.PrevAE[r,e,y] == model.AnnualEmissions[r,e,y]
        model.FixedAE_constraint = Constraint(model.REGION, model.EMISSION, model.FIXEDMILPYEARS, rule=FixedAE_rule)        
        
    def FixedConstraintCapacity_rule(model,r,l,t,m,y):
        return model.PrevROA[r,l,t,m,y] == model.RateOfActivity[r,l,t,m,y]
    model.FixedConstraintCapacity = Constraint(model.REGION, model.TIMESLICE, model.TECHNOLOGY, model.MODE_OF_OPERATION, model.FIXEDMILPYEARS, rule=FixedConstraintCapacity_rule)

    def FixedSalvageValueAtEndOfPeriod1_rule(model,r,t,y): 
        return model.PrevSV[r,t,y] == model.SalvageValue[r,t,y]       
    model.FixedSalvageValueAtEndOfPeriod1 = Constraint(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, rule=FixedSalvageValueAtEndOfPeriod1_rule)

    def FixedSalvageValueDiscountedToStartYear_rule(model,r,t,y):
        return model.PrevDSV[r,t,y] == model.DiscountedSalvageValue[r,t,y]
    model.FixedSalvageValueDiscountedToStartYear = Constraint(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, rule=FixedSalvageValueDiscountedToStartYear_rule)

    def FixedNewCap_rule(model,r,t,y): 
        return model.PrevNC[r,t,y] == model.NewCapacity[r,t,y]       
    model.FixedNewCap = Constraint(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, rule=FixedNewCap_rule)

    def FixedNNT_rule(model,r,t,y): 
        return model.PrevNNT[r,t,y] == model.NumberOfNewTechnologyUnits[r,t,y]     
    model.FixedNNT = Constraint(model.REGION, model.TECHNOLOGY, model.FIXEDMILPYEARS, rule=FixedNNT_rule)

    def FixedHy1_RateOfReservoirDischarge_rule(model,v,y,ls,r):
        return model.PrevRRD[v,y,ls,r] == model.RateOfReservoirDischarge[v,y,ls,r]
    model.FixedHy1_RateOfReservoirDischarge = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedHy1_RateOfReservoirDischarge_rule)
    
    def FixedHy9_and_Hy10_ReservoirLevelSeason_rule(model,v,y,ls,r): #controlla cosa sono quegli 1e-3 della riga 367-->è il /1000 che considera evap in mm/m2
        return model.PrevRLS[v,y,ls,r] == model.ReservoirLevelSeason[v,y,ls,r]
    model.FixedHy9_and_Hy10_ReservoirLevelSeason = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedHy9_and_Hy10_ReservoirLevelSeason_rule)

    def FixedSpillage_rule(model,v,y,ls,r):
        return model.PrevRRS[v,y,ls,r] == model.RateOfReservoirSpillage[v,y,ls,r]
    model.FixedSpillage = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedSpillage_rule)

    def FixedRelease_rule(model,v,y,ls,r):
        return model.PrevREL[v,y,ls,r] == model.Release[v,y,ls,r]
    model.FixedRelease = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedRelease_rule)

    def FixedInflow_rule(model,v,y,ls,r):
        return model.PrevINF[v,y,ls,r] == model.Inflow[v,y,ls,r]
    model.FixedInflow = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedInflow_rule)

    def FixedX_rule(model,v,y,ls,r):
        return model.PrevX[v,y,ls,r] == model.x[v,y,ls,r]
    model.FixedX = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedX_rule)

    def FixedY_rule(model,v,y,ls,r):
        return model.PrevY[v,y,ls,r] == model.y[v,y,ls,r]
    model.FixedY = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedY_rule)

    def FixedSurface_rule(model,v,y,ls,r):
        return model.PrevSUR[v,y,ls,r] == model.Surface[v,y,ls,r]
    model.FixedSurface = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedSurface_rule)

    def FixedReleaseLag_rule(model,v,y,ls,r):
        return model.PrevLAG[v,y,ls,r] == model.ReleaseLag[v,y,ls,r]
    model.FixedReleaseLag = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, model.REGION, rule=FixedReleaseLag_rule)

    def FixedMaxRel_rule(model,v,y,ls):
        return model.PrevMRE[v,y,ls] == model.MaxRel[v,y,ls]
    model.FixedMaxRel = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, rule=FixedMaxRel_rule)
    
    def FixedMinRel_rule(model,v,y,ls):
        return model.PrevmRE[v,y,ls] == model.MinRel[v,y,ls]
    model.FixedMinRel = Constraint(model.RESERVOIR, model.FIXEDMILPYEARS, model.SEASON, rule=FixedMinRel_rule)

###########################
#  Instantiate the model  #
###########################

instance = model.create_instance(filename=path[:-7]+"dat_files/FixAndRelax/"+PastOrFuture+"/"+nomefile+".dat",name="OS_FixAndRelax", report_timing=True)
instance.write(filename=path[:-7] +"lp_files/FixAndRelax/"+PastOrFuture+"/"+nomefile+".lp", io_options={'symbolic_solver_labels':True})

# SPECIFY YOUR CPLEX PATH
os.environ['PATH'] += os.pathsep + '/opt/ibm/ILOG/CPLEX_Studio128/cplex/bin/x86-64_linux'
import pyutilib
pyutilib.services.register_executable('cplex')

opt = SolverFactory('cplex')

option={}
for opzione in DATA.OpzioniCPLEX():
    option[opzione[0]]=opzione[1]
results = opt.solve(instance,tee=True, keepfiles=True,options=option)

#############################################################   

model.solutions.load_from(results)
with open(path[:-7] + "output_files/FixAndRelax/"+PastOrFuture+"/solution_"+nomefile+".txt","w") as f:
    instance.display(ostream=f)

end=datetime.datetime.now()
TELEGRAM.send([676406408],"SIMULINO OSE-FER MILPYEARS: " + str(MinYear)+"-"+str(MaxYear) +"   overlap: " +str(OverlappingSize)+ "   tempo: {}".format(end-start))

