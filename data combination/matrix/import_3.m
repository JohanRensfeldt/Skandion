function x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = import_3(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb)
	% Fill missing data
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    fillmissing(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","MaxGap",0.7,"DataVariables","Fx___");
	% Fill outliers
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    filloutliers(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "pchip","quartiles","ThresholdFactor",2.5,...
	    "DataVariables",["Fy___","Sx___","Sy___"]);
	% Fill missing data
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    fillmissing(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","MaxGap",10,"DataVariables",["x2DF___","x2DSL_R___","x2DSU_D___"]);
	% Fill outliers
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    filloutliers(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","movmedian",60,"ThresholdFactor",4.5,...
	    "DataVariables",["x2DF___","x2DSL_R___","x2DSU_D___"]);
	% Fill missing data
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    fillmissing(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","MaxGap",10,...
	    "DataVariables",["DCen_Gy_","PLx_mm_","PRx_mm_","PLy_mm_"]);
	% Fill outliers
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    filloutliers(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","movmedian",54,"ThresholdFactor",8,...
	    "DataVariables",["DCen_Gy_","PLx_mm_","PRx_mm_","PLy_mm_"]);
	% Fill missing data
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    fillmissing(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","MaxGap",10,...
	    "DataVariables",["PRy_mm_","FWHMX_mm_","FWHMY_mm_","Fx___1","Fy___1","Sx___1",...
	    "Sy___1","x2DF___1"]);
	% Fill outliers
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    filloutliers(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","movmedian",20,"ThresholdFactor",7,...
	    "DataVariables",["PRy_mm_","Fx___1","Fy___1","Sx___1","Sy___1","x2DF___1"]);
	% Fill missing data
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    fillmissing(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","MaxGap",10,...
	    "DataVariables",["x2DSL_R___1","x2DSU_D___1","PLx_mm_1","PRx_mm_1","PLy_mm_1",...
	    "PRy_mm_1"]);
	% Fill outliers
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    filloutliers(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","movmedian",33,"ThresholdFactor",7,...
	    "DataVariables",["x2DSL_R___1","x2DSU_D___1"]);
	% Fill outliers
	x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb = ...
	    filloutliers(x01_DailyCubeParametersMatriXXQA_referenceFromMeasurements_comb,...
	    "linear","movmedian",120,"ThresholdFactor",10,...
	    "DataVariables","DiffR15_SOBP10_2Gy");
end