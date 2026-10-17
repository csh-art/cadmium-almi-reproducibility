source(file.path("01_nhanes", "scripts", "00_config.R"), encoding = "UTF-8")
core_path <- file.path(generated_dir, "nhanes_survey_base_20_59_core.rds")
if (!file.exists(core_path)) stop("Run 01_build_core_dataset.R first.", call. = FALSE)
survey_base <- readRDS(core_path)

merge_component <- function(x, cyc, suf, stem, vars, special = NULL) {
  path <- make_path(cyc, stem, suf, special)
  d <- read_xpt_checked(path)
  for (v in vars) if (!v %in% names(d)) d[[v]] <- NA_real_
  d <- d[, unique(c("SEQN", vars)), drop = FALSE]
  merge(x, d, by = "SEQN", all.x = TRUE, sort = FALSE)
}

recode_alcohol <- function(x, cyc) {
  out <- rep(NA_character_, nrow(x))
  if (cyc != "2017-2018") {
    q <- x$ALQ120Q; u <- x$ALQ120U; screen <- x$ALQ101
    valid_q <- !is.na(q) & !q %in% c(777, 999)
    valid_u <- !is.na(u) & u %in% 1:3
    weekly <- valid_q & valid_u & (u == 1 | (u == 2 & q >= 5) | (u == 3 & q >= 52))
    less <- valid_q & valid_u & ((u == 2 & q >= 1 & q < 5) | (u == 3 & q >= 12 & q < 52))
    invalid_freq <- (!is.na(q) & q %in% c(777, 999)) | (!is.na(u) & !u %in% 1:3)
    missing <- is.na(screen) | invalid_freq
    out[!missing] <- "None or rare"
    out[less] <- "Less than weekly"
    out[weekly] <- "Weekly or more"
  } else {
    screen <- x$ALQ111; f <- x$ALQ121
    missing <- is.na(screen) | screen %in% c(7, 9) | (!is.na(f) & f %in% c(77, 99))
    out[!missing] <- "None or rare"
    out[!missing & f %in% c(6, 7)] <- "Less than weekly"
    out[!missing & f %in% 1:5] <- "Weekly or more"
  }
  factor(out, levels = c("None or rare", "Less than weekly", "Weekly or more"))
}

recode_activity <- function(x) {
  specs <- list(
    c("PAQ605", "PAQ610", "PAD615", 2),
    c("PAQ620", "PAQ625", "PAD630", 1),
    c("PAQ635", "PAQ640", "PAD645", 1),
    c("PAQ650", "PAQ655", "PAD660", 2),
    c("PAQ665", "PAQ670", "PAD675", 1)
  )
  total <- rep(0, nrow(x)); missing <- rep(FALSE, nrow(x))
  for (s in specs) {
    yn <- s[[1]]; days <- s[[2]]; mins <- s[[3]]; mult <- as.numeric(s[[4]])
    yes <- x[[yn]] == 1; no <- x[[yn]] == 2
    d <- x[[days]]; m <- x[[mins]]
    dvalid <- !is.na(d) & !d %in% c(77, 99, 777, 999)
    mvalid <- !is.na(m) & !m %in% c(777, 999, 7777, 9999)
    total[yes & dvalid & mvalid] <- total[yes & dvalid & mvalid] + d[yes & dvalid & mvalid] * m[yes & dvalid & mvalid] * mult
    missing <- missing | (!yes & !no) | (yes & !(dvalid & mvalid))
  }
  out <- ifelse(missing, NA_character_, ifelse(total >= 150, "Meets guideline", "Below guideline"))
  factor(out, levels = c("Meets guideline", "Below guideline"))
}

egfr_2021 <- function(scr, age, sex) {
  female <- sex == "Female"
  k <- ifelse(female, 0.7, 0.9)
  alpha <- ifelse(female, -0.241, -0.302)
  sex_factor <- ifelse(female, 1.012, 1)
  ratio <- scr / k
  142 * pmin(ratio, 1)^alpha * pmax(ratio, 1)^(-1.200) * 0.9938^age * sex_factor
}

out_parts <- list()
for (i in seq_len(nrow(cycle_spec))) {
  cyc <- cycle_spec$cycle[i]; suf <- cycle_spec$suffix[i]
  x <- survey_base[survey_base$survey_cycle == cyc, , drop = FALSE]
  x <- merge_component(x, cyc, suf, "SMQ", c("SMQ020", "SMQ040"))
  x <- merge_component(x, cyc, suf, "ALQ", c("ALQ101", "ALQ111", "ALQ120Q", "ALQ120U", "ALQ121", "ALQ130"))
  x <- merge_component(x, cyc, suf, "PAQ", c("PAQ605", "PAQ610", "PAD615", "PAQ620", "PAQ625", "PAD630", "PAQ635", "PAQ640", "PAD645", "PAQ650", "PAQ655", "PAD660", "PAQ665", "PAQ670", "PAD675"))
  x <- merge_component(x, cyc, suf, "BPQ", c("BPQ020", "BPQ080", "BPQ100D"))
  x <- merge_component(x, cyc, suf, "DIQ", c("DIQ010"))
  x <- merge_component(x, cyc, suf, "BIOPRO", c("LBXSCR"), special = paste0("BIOPRO_", suf, "_SUPP.XPT"))
  cot_file <- if (suf == "G") "COTNAL_G_SUPP.XPT" else paste0("COT_", suf, "_SUPP.XPT")
  x <- merge_component(x, cyc, suf, "COT", c("LBXCOT"), special = cot_file)

  x$age <- x$RIDAGEYR
  x$sex <- factor(ifelse(x$RIAGENDR == 1, "Male", ifelse(x$RIAGENDR == 2, "Female", NA_character_)), levels = c("Female", "Male"))
  race_map <- c(`1`="Mexican American", `2`="Other Hispanic", `3`="Non-Hispanic White", `4`="Non-Hispanic Black", `6`="Non-Hispanic Asian", `7`="Other or multiracial")
  x$race <- factor(unname(race_map[as.character(x$RIDRETH3)]), levels = c("Non-Hispanic White", "Non-Hispanic Black", "Mexican American", "Other Hispanic", "Non-Hispanic Asian", "Other or multiracial"))

  edu <- rep(NA_character_, nrow(x)); edu[x$DMDEDUC2 %in% c(1,2)] <- "Less than high school"; edu[x$DMDEDUC2 == 3] <- "High school or GED"; edu[x$DMDEDUC2 == 4] <- "Some college or associate degree"; edu[x$DMDEDUC2 == 5] <- "College graduate or above"
  x$education <- factor(edu, levels = c("College graduate or above", "Some college or associate degree", "High school or GED", "Less than high school"))
  mar <- rep(NA_character_, nrow(x)); mar[x$DMDMARTL %in% c(1,6)] <- "Married or living with partner"; mar[x$DMDMARTL == 5] <- "Never married"; mar[x$DMDMARTL %in% c(2,3,4)] <- "Widowed, divorced, or separated"
  x$marital <- factor(mar, levels = c("Married or living with partner", "Never married", "Widowed, divorced, or separated"))
  x$pir_continuous <- x$INDFMPIR

  sm <- rep(NA_character_, nrow(x)); sm[x$SMQ020 == 2] <- "Never"; sm[x$SMQ020 == 1 & x$SMQ040 == 3] <- "Former"; sm[x$SMQ020 == 1 & x$SMQ040 %in% c(1,2)] <- "Current"
  x$smoking_status <- factor(sm, levels = c("Never", "Former", "Current"))
  x$alcohol_status <- recode_alcohol(x, cyc)
  x$physical_activity <- recode_activity(x)
  x$hypertension <- factor(ifelse(x$BPQ020 == 1, "Yes", ifelse(x$BPQ020 == 2, "No", NA_character_)), levels = c("No", "Yes"))
  x$diabetes <- factor(ifelse(x$DIQ010 == 1, "Yes", ifelse(x$DIQ010 == 2, "No", ifelse(x$DIQ010 == 3, "Borderline", NA_character_))), levels = c("No", "Borderline", "Yes"))
  hyp_yes <- x$BPQ080 == 1 | x$BPQ100D == 1
  hyp_no <- !hyp_yes & x$BPQ080 == 2
  hyp <- rep(NA_character_, nrow(x)); hyp[hyp_yes] <- "Yes"; hyp[hyp_no] <- "No"
  x$hyperlipidemia <- factor(hyp, levels = c("No", "Yes"))

  x$serum_creatinine_calibrated <- x$LBXSCR
  if (cyc == "2017-2018") x$serum_creatinine_calibrated <- 1.051 * x$LBXSCR - 0.06945
  x$egfr_main <- egfr_2021(x$serum_creatinine_calibrated, x$age, x$sex)
  x$log2_serum_cotinine <- ifelse(!is.na(x$LBXCOT) & x$LBXCOT > 0, log2(x$LBXCOT), NA_real_)
  x$bmi <- x$BMXBMI
  x$waist_circumference_cm <- x$BMXWAIST
  out_parts[[i]] <- x
}

final_survey_base <- do.call(rbind, out_parts)
strict <- final_survey_base[final_survey_base$core_analysis_eligible_strict, , drop = FALSE]

expected_nonmissing <- c(age=7628, sex=7628, race=7628, education=7627, marital=7627, pir_continuous=7012,
                         smoking_status=7624, alcohol_status=7170, physical_activity=7588, hypertension=7621,
                         diabetes=7622, hyperlipidemia=7594, egfr_main=7497, log2_serum_cotinine=7539,
                         log2_blood_lead=7628, bmi=7624, waist_circumference_cm=7577)
observed <- vapply(names(expected_nonmissing), function(v) sum(!is.na(strict[[v]])), integer(1))
if (!all(observed == expected_nonmissing)) {
  stop("Covariate nonmissing-count validation failed:\n", paste(names(observed), observed, expected_nonmissing, sep = ": ", collapse = "\n"), call. = FALSE)
}

saveRDS(final_survey_base, file.path(generated_dir, "nhanes_survey_base_20_59_final.rds"), compress = "xz")
saveRDS(strict, file.path(generated_dir, "nhanes_analysis_dataset_20_59_pre_model.rds"), compress = "xz")
utils::write.csv(data.frame(variable=names(observed), nonmissing_n=observed, expected_n=expected_nonmissing), file.path(generated_dir, "covariate_nonmissing_validation.csv"), row.names = FALSE)
cat("02_build_covariates.R PASS\n")
