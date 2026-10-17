source(file.path("01_nhanes", "scripts", "00_config.R"), encoding = "UTF-8")
p <- file.path(generated_dir,"survey_base_with_quartile_trend_score.rds")
if (!file.exists(p)) stop("Run 04_weighted_quartiles.R first.",call.=FALSE)
d <- readRDS(p)
design <- survey::svydesign(ids=~survey_psu_id,strata=~survey_strata_id,weights=~weight_8yr,data=d,nest=TRUE)

m3_cov <- c("age","sex","race","survey_cycle","education","marital","pir_continuous","smoking_status","alcohol_status","physical_activity")
m4_cov <- c(m3_cov,"hypertension","diabetes","hyperlipidemia")
fit_one <- function(des,covars,drop_covars=character()) {
 covars <- setdiff(covars,drop_covars); f <- as.formula(paste("ALMI_strict ~ log2_blood_cadmium",if(length(covars)) paste("+",paste(covars,collapse=" + ")) else "")); fit <- survey::svyglm(f,design=des); co <- summary(fit)$coefficients; est <- co["log2_blood_cadmium","Estimate"]; se <- co["log2_blood_cadmium","Std. Error"]; ddf <- survey::degf(des); crit <- qt(.975,ddf); data.frame(estimate=est,se=se,lower_95ci=est-crit*se,upper_95ci=est+crit*se,p_value=co["log2_blood_cadmium",ncol(co)],design_df=ddf)
}
rows <- list()
# Additional adjustment, based on Model 4
for (nm in c("S1_eGFR","S2_Cotinine","S3_Lead","S4_Full")) {
 flag <- paste0("eligible_",nm); extra <- switch(nm,S1_eGFR="egfr_main",S2_Cotinine="log2_serum_cotinine",S3_Lead="log2_blood_lead",S4_Full=c("egfr_main","log2_serum_cotinine","log2_blood_lead")); des <- subset(design,get(flag)); rows[[nm]] <- cbind(analysis=nm,N=sum(d[[flag]]),fit_one(des,c(m4_cov,extra)))
}
# Common Model 4 sample
common <- d$eligible_M4_Clinical
des_common <- subset(design,eligible_M4_Clinical)
for (nm in c("M1","M2","M3","M4")) {
 cv <- switch(nm,M1=character(),M2=c("age","sex","race","survey_cycle"),M3=m3_cov,M4=m4_cov); rows[[paste0("Common_",nm)]] <- cbind(analysis=paste0("Common_",nm),N=sum(common),fit_one(des_common,cv))
}
# Restricted samples based on M3 eligibility
p99 <- as.numeric(survey::svyquantile(~LBXBCD,subset(design,eligible_M3_Primary),.99,ci=FALSE,na.rm=TRUE))
rows[["Exclude_current"]] <- cbind(analysis="Exclude_current_smokers",N=sum(d$eligible_M3_Primary & d$smoking_status!="Current"),fit_one(subset(design,eligible_M3_Primary & smoking_status!="Current"),m3_cov))
rows[["Never"]] <- cbind(analysis="Never_smokers_only",N=sum(d$eligible_M3_Primary & d$smoking_status=="Never"),fit_one(subset(design,eligible_M3_Primary & smoking_status=="Never"),m3_cov,drop_covars="smoking_status"))
rows[["P99"]] <- cbind(analysis="Exclude_highest_1pct_cadmium",N=sum(d$eligible_M3_Primary & d$LBXBCD<=p99),fit_one(subset(design,eligible_M3_Primary & LBXBCD<=p99),m3_cov))
# Sex specific and interaction
for (sx in levels(d$sex)) rows[[paste0("Sex_",sx)]] <- cbind(analysis=paste0("Sex_",sx),N=sum(d$eligible_M3_Primary & d$sex==sx),fit_one(subset(design,eligible_M3_Primary & sex==sx),m3_cov,drop_covars="sex"))
des_m3 <- subset(design,eligible_M3_Primary); intfit <- survey::svyglm(as.formula(paste("ALMI_strict ~ log2_blood_cadmium * sex +",paste(setdiff(m3_cov,"sex"),collapse=" + "))),design=des_m3); ip <- survey::regTermTest(intfit,~log2_blood_cadmium:sex)$p
# Leave one cycle out
for(cyc in levels(factor(d$survey_cycle))) rows[[paste0("LOCO_",cyc)]] <- cbind(analysis=paste0("Exclude_",cyc),N=sum(d$eligible_M3_Primary & d$survey_cycle!=cyc),fit_one(subset(design,eligible_M3_Primary & survey_cycle!=cyc),m3_cov))
res <- do.call(rbind,rows); rownames(res)<-NULL
utils::write.csv(res,file.path(generated_dir,"sensitivity_effect_results.csv"),row.names=FALSE)
utils::write.csv(data.frame(interaction="blood_cadmium_by_sex",p_value=ip),file.path(generated_dir,"sex_interaction_result.csv"),row.names=FALSE)
utils::write.csv(data.frame(weighted_p99_cadmium=p99),file.path(generated_dir,"cadmium_p99.csv"),row.names=FALSE)
if (abs(ip-0.731948979532676)>5e-6) warning("Sex-interaction p differs slightly from archived value; inspect factor coding/software version.")
cat("06_sensitivity_analyses.R completed.\n")
