
**** Statistical analysis code for intention-to-treat analysis

// Pooled logistic regression model was fitted for each outcome:

foreach outcome in cvd mi hf stroke death muscleAE liver{

	logit `outcome' statin_initiate /*
		*/i.sex /*
		*/age_bl /*
		*/hba1c_bl /*
		*/sbp_bl /*
		*/dbp_bl /*
		*/ldl_c_bl /*
		*/hdl_c_bl /*
		*/tc_bl /*
		*/egfr_bl /*	
		*/cci_bl /*
		*/i.hypertension_bl /*
		*/i.obesity_bl /*
		*/i.pvd_bl /*
		*/i.af_bl /*
		*/i.copd_bl /*
		*/i.renal_bl /*
		*/i.dementia_bl /*
		*/i.acei_arb_bl /*
		*/i.b_blocker_bl /*
		*/i.ccb_bl /*
		*/i.diuretic_bl /*
		*/i.oral_dm_drug_bl /*
		*/i.insulin_bl /*
		*/i.aspirin_bl /*
		*/i.sopc_attend_bl /*
		*/i.hosp_bl /*
		*/i.smoking_bl /*
		*/time /*
		*/time_square /*
		*/i.baseline  if time<=fup_`outcome'_obs, cluster(id_seq) or

}
