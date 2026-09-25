*ModelSim.do
*-------------------v.06.12.2019; @A. Michaud for AddedWorker w/ K.Ellieroth---------------------------------*

clear all
set more off

*Set paths
foreach v in "Baseline" "Cohort1950" "Cohort1960"  "Cohort1970"  "Cohort1980"   "Kap_dcr" "RoE_incr" "Wgap_dcr"  {
global home_dir "C:\Users\IRAMM03\Desktop\KEAM\Matlab"
global data_dir "$home_dir\Output"
global simdata_dir "`v'"
global stata_dir "$home_dir\STATA"
global simcsv_dir "SimPanel"
global figOUT_dir "Figures"
global tabOUT_dir "Tables"

*-------------------------------------------------------------------------------------------------------*
*This file analyzes model simulated data and compares to empirical data from the CPS and PSID
*--------------------------------------------------------------------------------------------------------

cd $data_dir
cd "$simdata_dir"

clear
import delimited "simdata.csv"

gen version="`v'"

rename v1 id
rename v2 time
rename v3 year
rename v4 quarter
rename v5 recession
rename v6 age
rename v7 employed
rename v8 hours
rename v9 wageAll
	label var wageAll "Wages including shadow wage"
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
rename v24 ValueFn
rename v25 UtilC
rename v26 UtilL

cd $stata_dir 
cd $simcsv_dir

save "`v'_sim" , replace
}
*-------------------------------------------
*Merge Data
cd $stata_dir 
cd $simcsv_dir
use "Baseline_sim", clear

foreach v in  "Cohort1950_sim" "Cohort1960_sim"  "Cohort1970_sim"  "Cohort1980_sim"  "Kap_dcr_sim" "RoE_incr_sim" "Wgap_dcr_sim"  {
	append using "`v'.dta"
}

save "Merged_sim" , replace
*------------------------------------------------------------------------

use "Merged_sim", clear
*------------------------------------------------------------------------
*its a panel
	egen panelid = group(id version), label
	xtset panelid time
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
	by panelid: egen ftkappa = max(ftKappa1)
	label var ftkappa "Individual fixed cost of work"
	drop ftKappa1
	
*Variable cost of work is when young
	gen Kap_y1 = workcost-ftkappa if age<40 
	by panelid: egen Kap_y = max(Kap_y1)
	label var Kap_y "Individual cost draw when young"	
	drop Kap_y1
	
replace recession=0 if recession==1
replace recession=1 if recession==2

gen ageD1=(age>25 & age<40)
gen ageD2=(age>39 & age<55)
gen ageD3=(age>54 & age<65)

gen yob= year-age
gen cohort=1 if yob>1933 & yob<1939
	replace cohort=2 if yob>1943 & yob<1949
	replace cohort=3 if yob>1953 & yob<1959
	replace cohort=4 if yob>1963 & yob<1969
	replace cohort=5 if yob>1973 & yob<1979	
	
gen PTemp=(employed==1 & FTemp~=1)	


*Calculate comp stats:
 *--Individual level:	
	*Income risk
	 by panelid: egen Incvar=sd(hh_inc)
	 by panelid: egen Expinc=mean(hh_inc) if recession==0
	 by panelid: egen Recinc=mean(hh_inc) if recession==1
	 by panelid: egen Hunempinc=mean(hh_inc) if h_emp==0
	 by panelid: egen Hempinc=mean(hh_inc) if h_emp==1
	*Utility 
		gen kapval = 0 if employed==0
			replace kapval = workcost if employed==1

	preserve 
		collapse (max) Incvar Expinc Recinc Hunempinc Hempinc (mean) UtilC UtilL kapval ValueFn, by(id version)
			gen d_bcinc = (Recinc-Expinc)/Expinc
			gen d_uinc = (Hunempinc-Hempinc)/Hempinc
		*Compare to baseline	
		 gen base=(version=="Baseline")
		 sort id
			foreach var of varlist Incvar Expinc Recinc Hunempinc Hempinc d_bcinc d_uinc UtilC UtilL kapval ValueFn {
				gen b = `var' if base==1
				by id: egen base_`var' = max(b)
				gen vdif_`var' = (`var'-base_`var')/base_`var'
				gen p90_`var' = vdif_`var'
				gen p10_`var' = vdif_`var'
				drop b
			}
		
		collapse (mean) base_* vdif_* (p90) p90_* (p10) p10_*, by(version)
			export delimited using "Individual_relto_Base", replace
	restore
	
 *--Aggregate Level
	gen employment = employed + h_emp
	gen hh_hours= hours + h_emp*40
	
	collapse (sum) employment hh_hours hh_inc, by(recession version)
	
	 gen base=(version=="Baseline")
	 
 			foreach var of varlist employment hh_hours hh_inc {
				gen rr=`var' if recession==0
				egen rec = max(rr)
				gen dbc_`var' = `var'-rec
				gen b = `var' if base==1
				gen bb = dbc_`var' if base==1
				egen base_dbc_`var' = max(bb)
				egen base_`var' = max(b)
				gen vdifbc_`var' = dbc_`var'-base_dbc_`var'
				gen vdif_`var' = (`var'-base_`var')/base_`var'
				drop b bb rec rr
			}
			
			keep if recession==1
			keep version vdif* base_*
			
	export delimited using "Agg_relto_Base", replace
 
 
	