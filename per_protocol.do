
**** Statistical analysis code for per-protocol analysis

//  Weighting procedure and the model fitting for per-protocol effect estimate were carried out for each outcome:

** 1. For the outcome of all-cause mortality

foreach outcome in death{
	
	* Model to predict the probablity of receiving statin therapy in treatment arm 	
	
	logit statin_using i.sex age_bl hba1c_bl sbp_bl dbp_bl ldl_c_bl hdl_c_bl tc_bl egfr_bl cci_bl i.hypertension_bl i.obesity_bl i.pvd_bl i.af_bl i.copd_bl i.renal_bl i.dementia_bl i.acei_arb_bl i.b_blocker_bl i.ccb_bl i.diuretic_bl i.oral_dm_drug_bl i.insulin_bl i.aspirin_bl i.sopc_attend_bl i.hosp_bl i.smoking_bl age hba1c sbp dbp ldl_c hdl_c tc egfr cci i.hypertension i.obesity i.pvd i.af i.copd i.renal i.dementia i.acei_arb i.b_blocker i.ccb i.diuretic i.insulin i.aspirin i.oral_dm_drug i.sopc_attend i.hosp i.smoke time time_square if model_`outcome'==1, cluster(id_seq)
	predict pd_statin_`outcome'_1 if model_`outcome'==1, pr

	* Model to predict the probablity of receiving statin therapy in control arm 	
	
	logit statin_using i.sex age_bl hba1c_bl sbp_bl dbp_bl ldl_c_bl hdl_c_bl tc_bl egfr_bl cci_bl i.hypertension_bl i.obesity_bl i.pvd_bl i.af_bl i.copd_bl i.renal_bl i.dementia_bl i.acei_arb_bl i.b_blocker_bl i.ccb_bl i.diuretic_bl i.oral_dm_drug_bl i.insulin_bl i.aspirin_bl i.sopc_attend_bl i.hosp_bl i.smoking_bl age hba1c sbp dbp ldl_c hdl_c tc egfr cci i.hypertension i.obesity i.pvd i.af i.copd i.renal i.dementia i.acei_arb i.b_blocker i.ccb i.diuretic i.insulin i.aspirin i.oral_dm_drug i.sopc_attend i.hosp i.smoke time time_square if model_`outcome'==2, cluster(id_seq)
	predict pd_statin_`outcome'_2 if model_`outcome'==2, pr

	* construction of IP weights for treatment adherence
	
	gen pd_statin_`outcome'= .
	replace pd_statin_`outcome'=pd_statin_`outcome'_1 if model_`outcome'==1
	replace pd_statin_`outcome'=pd_statin_`outcome'_2 if model_`outcome'==2

	gen weight_point_`outcome'= .
	replace weight_point_`outcome'=1/pd_statin_`outcome' if model_`outcome'==1
	replace weight_point_`outcome'=1/(1-pd_statin_`outcome') if model_`outcome'==2
	replace weight_point_`outcome'=1 if time>fup_`outcome'_weight
	replace weight_point_`outcome'=1 if time==0
	
	sort id_seq time
	gen weight_cum_`outcome'=weight_point_`outcome'
	by id_seq: replace weight_cum_`outcome'=weight_cum_`outcome'*weight_cum_`outcome'[_n-1] if _n!=1

	// truncate the outliers of IP weight at 10
	gen sw_`outcome'=weight_cum_`outcome'
	replace sw_`outcome'=10 if sw_`outcome'>10

	* Model for per-protocol effect estimate
	
	logit `outcome' arm i.sex age_bl hba1c_bl sbp_bl dbp_bl ldl_c_bl hdl_c_bl tc_bl egfr_bl cci_bl i.hypertension_bl i.obesity_bl i.pvd_bl i.af_bl i.copd_bl i.renal_bl i.dementia_bl i.acei_arb_bl i.b_blocker_bl i.ccb_bl i.diuretic_bl i.oral_dm_drug_bl i.insulin_bl i.aspirin_bl i.sopc_attend_bl i.hosp_bl i.smoking_bl time time_square i.baseline [pw=sw_`outcome'] if time<=fup_`outcome', cluster(id_seq) or

}

** 2. For other outcomes: additionally adjusted by a time-varying weight of not dying to account for the competing risk

foreach outcome in cvd mi hf stroke liver muscleAE{

	* Model to predict the probablity of receiving statin therapy in treatment arm 	

	logit statin_using i.sex age_bl hba1c_bl sbp_bl dbp_bl ldl_c_bl hdl_c_bl tc_bl egfr_bl cci_bl i.hypertension_bl i.obesity_bl i.pvd_bl i.af_bl i.copd_bl i.renal_bl i.dementia_bl i.acei_arb_bl i.b_blocker_bl i.ccb_bl i.diuretic_bl i.oral_dm_drug_bl i.insulin_bl i.aspirin_bl i.sopc_attend_bl i.hosp_bl i.smoking_bl age hba1c sbp dbp ldl_c hdl_c tc egfr cci i.hypertension i.obesity i.pvd i.af i.copd i.renal i.dementia i.acei_arb i.b_blocker i.ccb i.diuretic i.insulin i.aspirin i.oral_dm_drug i.sopc_attend i.hosp i.smoke time time_square if model_`outcome'==1, cluster(id_seq)
	predict pd_statin_`outcome'_1 if model_`outcome'==1, pr

	* Model to predict the probablity of receiving statin therapy in control arm 	

	logit statin_using i.sex age_bl hba1c_bl sbp_bl dbp_bl ldl_c_bl hdl_c_bl tc_bl egfr_bl cci_bl i.hypertension_bl i.obesity_bl i.pvd_bl i.af_bl i.copd_bl i.renal_bl i.dementia_bl i.acei_arb_bl i.b_blocker_bl i.ccb_bl i.diuretic_bl i.oral_dm_drug_bl i.insulin_bl i.aspirin_bl i.sopc_attend_bl i.hosp_bl i.smoking_bl age hba1c sbp dbp ldl_c hdl_c tc egfr cci i.hypertension i.obesity i.pvd i.af i.copd i.renal i.dementia i.acei_arb i.b_blocker i.ccb i.diuretic i.insulin i.aspirin i.oral_dm_drug i.sopc_attend i.hosp i.smoke time time_square if model_`outcome'==2, cluster(id_seq)
	predict pd_statin_`outcome'_2 if model_`outcome'==2, pr

	*** Model to predict the probablity of death
	
	logit death arm i.sex age_bl hba1c_bl sbp_bl dbp_bl ldl_c_bl hdl_c_bl tc_bl egfr_bl cci_bl i.hypertension_bl i.obesity_bl i.pvd_bl i.af_bl i.copd_bl i.renal_bl i.dementia_bl i.acei_arb_bl i.b_blocker_bl i.ccb_bl i.diuretic_bl i.oral_dm_drug_bl i.insulin_bl i.aspirin_bl i.sopc_attend_bl i.hosp_bl i.smoking_bl age hba1c sbp dbp ldl_c hdl_c tc egfr cci i.hypertension i.obesity i.pvd i.af i.copd i.renal i.dementia i.acei_arb i.b_blocker i.ccb i.diuretic i.insulin i.aspirin i.oral_dm_drug i.sopc_attend i.hosp i.smoke time time_square if model_compete==1, cluster(id_seq)
	predict pd_death_`outcome' if model_compete==1, pr
	
	* construction of IP weights for treatment adherence and remaining alive

	gen pd_statin_`outcome'= .
	replace pd_statin_`outcome'=pd_statin_`outcome'_1 if model_`outcome'==1
	replace pd_statin_`outcome'=pd_statin_`outcome'_2 if model_`outcome'==2

	gen weight_point_`outcome'= .
	replace weight_point_`outcome'=1/pd_statin_`outcome' if model_`outcome'==1
	replace weight_point_`outcome'=1/(1-pd_statin_`outcome') if model_`outcome'==2
	replace weight_point_`outcome'=1 if time>fup_`outcome'_weight

	gen weight_nodying_`outcome'=1/(1-pd_death_`outcome') if model_compete==1
	
	gen weight_point_`outcome'_adj = weight_point_`outcome'*weight_nodying_`outcome'
	replace weight_point_`outcome'_adj=1 if time==0

	sort id_seq time
	gen weight_cum_`outcome'=weight_point_`outcome'_adj
	by id_seq: replace weight_cum_`outcome'=weight_cum_`outcome'*weight_cum_`outcome'[_n-1] if _n!=1

	// truncate the outliers of IP weight at 10
	gen sw_`outcome'=weight_cum_`outcome'
	replace sw_`outcome'=10 if sw_`outcome'>10
	
	* Model for per-protocol effect estimate
	
	logit outcome_`outcome' arm i.sex age_bl hba1c_bl sbp_bl dbp_bl ldl_c_bl hdl_c_bl tc_bl egfr_bl cci_bl i.hypertension_bl i.obesity_bl i.pvd_bl i.af_bl i.copd_bl i.renal_bl i.dementia_bl i.acei_arb_bl i.b_blocker_bl i.ccb_bl i.diuretic_bl i.oral_dm_drug_bl i.insulin_bl i.aspirin_bl i.sopc_attend_bl i.hosp_bl i.smoking_bl time time_square i.baseline [pw=sw_`outcome'] if time<=fup_`outcome', cluster(id_seq) or


}






