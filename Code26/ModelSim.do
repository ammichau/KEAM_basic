*ModelSim.do
*-------------------v.06.12.2019; @A. Michaud for AddedWorker w/ K.Ellieroth---------------------------------*

clear all
set more off

*Set paths
cd "C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab"

global HOME_dir "C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab"
global rawdata_dir "C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab\Output"
global figOUT_dir "Figures"
global tabOUT_dir "Tables"

*-------------------------------------------------------------------------------------------------------*
*This file analyzes model simulated data and compares to empirical data from the CPS and PSID
*--------------------------------------------------------------------------------------------------------

import delimited C:\Users\amichau9\Dropbox\DemogBCtrendsCode\Code_Simple_Model\Matlab\Output\simdata.csv

rename v1 id
rename v2 time
rename v3 year
rename v4 quarter
rename v5 recession
rename v6 age
rename v7 employed
rename v8 hours
rename v9 wageAll
	label var searchin "Wages including shadow wage"
rename v10 searchin
	label var searchin "Search Intensity"
rename v11 experience
rename v12 workcost
	label var workcost "Fixed Cost of Work (kappa)"
rename v13 ftWage
	label var ftWage "Wage fixed type"
rename v14 hh_inc
	label var hh_inc "Household Income"
rename v15 h_stat
	label var h_stat "Husband Status"
	label define h_stat 1 "employed" 2 "recentU" 3 "unemployed"
rename v16 jobloss
	label var jobloss "Involuntary Separation"
rename v17 quit
	label var quit "Voluntary Quit"
rename v18 FTemp
rename v19 nilf
rename v20 PTwoman
rename v21 cyclewoman
rename v22 careerwoman
rename v23 NiLFwoman

*------------------------------------------------------------------------
*its a panel
	xtset id time
*------------------------------------------------------------------------
*Be careful with these variables, need to coordinate for now
	
*------------------------------------------------------------------------
*Set up remaining variables

gen wage=wageAll if employed==1

gen h_emp = 1 if h_stat==1 | h_stat==2
	replace h_emp = 0 if h_stat==3
	label var h_emp "Husband Employed"
	
*Only have fixed cost of work when 40+	
gen ftKappa1 = workcost if age> 39 
	by id: egen ftkappa = max(ftKappa1)
	label var ftkappa "Individual fixed cost of work"
	drop ftKappa1
	
*Variable cost of work is when young
	gen Kap_y1 = workcost-ftkappa if age<40 
	by id: egen Kap_y = max(Kap_y1)
	label var Kap_y "Individual cost draw when young"	
	drop Kap_y1
	
replace recession=0 if recession==1
replace recession=1 if recession==2

gen ageD1=(Age>25 & Age<40)
gen ageD2=(Age>39 & Age<55)
gen ageD3=(Age>54 & Age<65)

gen yob= year-age
gen cohort=1 if yob>1933 & yob<1939
	replace cohort=2 if yob>1943 & yob<1949
	replace cohort=3 if yob>1953 & yob<1959
	replace cohort=4 if yob>1963 & yob<1969
	replace cohort=5 if yob>1973 & yob<1979	
	
gen PTemp=(employed==1 & FTemp~=1)	
	
	
	