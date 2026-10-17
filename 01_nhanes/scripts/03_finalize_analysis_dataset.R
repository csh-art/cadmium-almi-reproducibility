source(file.path("01_nhanes", "scripts", "00_config.R"), encoding = "UTF-8")
path <- file.path(generated_dir, "nhanes_survey_base_20_59_final.rds")
if (!file.exists(path)) stop("Run 02_build_covariates.R first.", call. = FALSE)
d <- readRDS(path)

m2 <- c("age","sex","race","survey_cycle")
m3 <- c(m2,"education","marital","pir_continuous","smoking_status","alcohol_status","physical_activity")
m4 <- c(m3,"hypertension","diabetes","hyperlipidemia")
s1 <- c(m4,"egfr_main")
s2 <- c(m4,"log2_serum_cotinine")
s3 <- c(m4,"log2_blood_lead")
s4 <- c(m4,"egfr_main","log2_serum_cotinine","log2_blood_lead")

cc <- function(vars) complete.cases(d[, vars, drop = FALSE])
d$eligible_M1_Crude <- d$core_analysis_eligible_strict
d$eligible_M2_Demographic <- d$core_analysis_eligible_strict & cc(m2)
d$eligible_M3_Primary <- d$core_analysis_eligible_strict & cc(m3)
d$eligible_M4_Clinical <- d$core_analysis_eligible_strict & cc(m4)
d$eligible_S1_eGFR <- d$core_analysis_eligible_strict & cc(s1)
d$eligible_S2_Cotinine <- d$core_analysis_eligible_strict & cc(s2)
d$eligible_S3_Lead <- d$core_analysis_eligible_strict & cc(s3)
d$eligible_S4_Full <- d$core_analysis_eligible_strict & cc(s4)

expected <- c(M1_Crude=7628, M2_Demographic=7628, M3_Primary=6575, M4_Clinical=6536, S1_eGFR=6423, S2_Cotinine=6460, S3_Lead=6536, S4_Full=6420)
flags <- paste0("eligible_", names(expected))
observed <- vapply(flags, function(v) sum(d[[v]]), integer(1)); names(observed) <- names(expected)
if (!all(observed == expected)) stop("Model sample-size validation failed.", call. = FALSE)

plan <- data.frame(
  model=c("M1_Crude","M2_Demographic","M3_Primary","M4_Clinical","S1_eGFR","S2_Cotinine","S3_Lead","S4_Full"),
  display_name=c("Model 1","Model 2","Model 3","Model 4","Sensitivity 1","Sensitivity 2","Sensitivity 3","Sensitivity 4"),
  covariates=c("",paste(m2,collapse=" + "),paste(m3,collapse=" + "),paste(m4,collapse=" + "),paste(s1,collapse=" + "),paste(s2,collapse=" + "),paste(s3,collapse=" + "),paste(s4,collapse=" + ")),
  stringsAsFactors=FALSE
)

saveRDS(d, file.path(generated_dir, "nhanes_survey_base_20_59_final_with_flags.rds"), compress="xz")
saveRDS(d[d$core_analysis_eligible_strict, , drop=FALSE], file.path(generated_dir, "nhanes_analysis_dataset_20_59_final.rds"), compress="xz")
utils::write.csv(plan, file.path(generated_dir, "fixed_analysis_model_plan.csv"), row.names=FALSE)
utils::write.csv(data.frame(model=names(observed), complete_case_n=as.integer(observed), expected_n=as.integer(expected)), file.path(generated_dir,"model_sample_size_validation.csv"), row.names=FALSE)
cat("03_finalize_analysis_dataset.R PASS\n")
