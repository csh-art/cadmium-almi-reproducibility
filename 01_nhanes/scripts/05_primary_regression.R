source(file.path("01_nhanes", "scripts", "00_config.R"), encoding = "UTF-8")
p <- file.path(generated_dir,"survey_base_with_quartile_trend_score.rds")
if (!file.exists(p)) stop("Run 04_weighted_quartiles.R first.", call.=FALSE)
d <- readRDS(p)
design <- survey::svydesign(ids=~survey_psu_id, strata=~survey_strata_id, weights=~weight_8yr, data=d, nest=TRUE)

models <- list(
 M1_Crude=list(flag="eligible_M1_Crude",covars=character()),
 M2_Demographic=list(flag="eligible_M2_Demographic",covars=c("age","sex","race","survey_cycle")),
 M3_Primary=list(flag="eligible_M3_Primary",covars=c("age","sex","race","survey_cycle","education","marital","pir_continuous","smoking_status","alcohol_status","physical_activity")),
 M4_Clinical=list(flag="eligible_M4_Clinical",covars=c("age","sex","race","survey_cycle","education","marital","pir_continuous","smoking_status","alcohol_status","physical_activity","hypertension","diabetes","hyperlipidemia"))
)

extract_term <- function(fit, term, design_df) {
  co <- summary(fit)$coefficients
  est <- unname(co[term,"Estimate"]); se <- unname(co[term,"Std. Error"]); p <- unname(co[term, ncol(co)])
  crit <- qt(.975, df=design_df)
  c(estimate=est,se=se,lower_95ci=est-crit*se,upper_95ci=est+crit*se,p_value=p)
}

cont_rows <- list(); quart_rows <- list(); test_rows <- list(); fit_store <- list()
for (nm in names(models)) {
  sp <- models[[nm]]; des <- subset(design, get(sp$flag)); ddf <- survey::degf(des)
  rhs <- c("log2_blood_cadmium",sp$covars); f <- as.formula(paste("ALMI_strict ~",paste(rhs,collapse=" + ")))
  fit <- survey::svyglm(f,design=des); e <- extract_term(fit,"log2_blood_cadmium",ddf)
  cont_rows[[nm]] <- data.frame(model=nm,N=sum(d[[sp$flag]]),design_df=ddf,t(e),row.names=NULL)
  fq <- as.formula(paste("ALMI_strict ~ cd_quartile +",paste(sp$covars,collapse=" + ")))
  if (!length(sp$covars)) fq <- ALMI_strict ~ cd_quartile
  qfit <- survey::svyglm(fq,design=des)
  qco <- summary(qfit)$coefficients
  for (term in c("cd_quartileQ2","cd_quartileQ3","cd_quartileQ4")) {
    e2 <- extract_term(qfit,term,ddf); quart_rows[[paste(nm,term)]] <- data.frame(model=nm,contrast=term,N=sum(d[[sp$flag]]),design_df=ddf,t(e2),row.names=NULL)
  }
  ft <- as.formula(paste("ALMI_strict ~ cd_quartile_median_score +",paste(sp$covars,collapse=" + ")))
  if (!length(sp$covars)) ft <- ALMI_strict ~ cd_quartile_median_score
  tfit <- survey::svyglm(ft,design=des); trend <- extract_term(tfit,"cd_quartile_median_score",ddf)
  overall <- survey::regTermTest(qfit,~cd_quartile)
  test_rows[[nm]] <- data.frame(model=nm,overall_quartile_p=overall$p,trend_p=trend[["p_value"]])
  fit_store[[nm]] <- list(continuous=fit,quartile=qfit,trend=tfit)
}
cont <- do.call(rbind,cont_rows); quart <- do.call(rbind,quart_rows); tests <- do.call(rbind,test_rows)
utils::write.csv(cont,file.path(generated_dir,"continuous_cadmium_results_all_models.csv"),row.names=FALSE)
utils::write.csv(quart,file.path(generated_dir,"quartile_results_long.csv"),row.names=FALSE)
utils::write.csv(tests,file.path(generated_dir,"quartile_overall_and_trend_tests.csv"),row.names=FALSE)
saveRDS(fit_store,file.path(generated_dir,"primary_regression_models.rds"),compress="xz")

m3 <- cont[cont$model=="M3_Primary",]
if (abs(m3$estimate - (-0.159643851702317)) > 5e-8) stop("Primary beta validation failed.",call.=FALSE)
cat("05_primary_regression.R PASS\n")
