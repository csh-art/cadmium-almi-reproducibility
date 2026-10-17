source(file.path("01_nhanes", "scripts", "00_config.R"), encoding = "UTF-8")

build_cycle_core <- function(spec) {
  cyc <- spec$cycle; suf <- spec$suffix; wvar <- spec$weight_variable
  demo <- read_xpt_checked(make_path(cyc, "DEMO", suf))
  pbcd <- read_xpt_checked(make_path(cyc, "PBCD", suf))
  dxx  <- read_xpt_checked(make_path(cyc, "DXX",  suf))
  bmx  <- read_xpt_checked(make_path(cyc, "BMX",  suf))

  demo$survey_cycle <- cyc
  demo$pbcd_available <- demo$SEQN %in% pbcd$SEQN

  keep_pbcd <- intersect(c("SEQN", "LBXBCD", "LBDBCDLC", "LBXBPB", "LBDBPBLC", "WTSH2YR"), names(pbcd))
  keep_dxx <- c("SEQN", "DXAEXSTS", "DXALATV", "DXARATV", "DXALLTV", "DXARLTV",
                "DXDLALE", "DXDRALE", "DXDLLLE", "DXDRLLE")
  keep_bmx <- intersect(c("SEQN", "BMXHT", "BMXWT", "BMXBMI", "BMXWAIST"), names(bmx))

  x <- merge(demo, pbcd[, keep_pbcd, drop = FALSE], by = "SEQN", all.x = TRUE, sort = FALSE)
  x <- merge(x, dxx[, keep_dxx, drop = FALSE], by = "SEQN", all.x = TRUE, sort = FALSE)
  x <- merge(x, bmx[, keep_bmx, drop = FALSE], by = "SEQN", all.x = TRUE, sort = FALSE)

  x$weight_2yr_raw <- if (wvar == "WTMEC2YR") x$WTMEC2YR else x$WTSH2YR
  x$weight_8yr <- x$weight_2yr_raw / 4
  x$log2_blood_cadmium <- ifelse(!is.na(x$LBXBCD) & x$LBXBCD > 0, log2(x$LBXBCD), NA_real_)
  x$log2_blood_lead <- ifelse(!is.na(x$LBXBPB) & x$LBXBPB > 0, log2(x$LBXBPB), NA_real_)

  tissue_vars <- c("DXALATV", "DXARATV", "DXALLTV", "DXARLTV")
  lean_vars <- c("DXDLALE", "DXDRALE", "DXDLLLE", "DXDRLLE")
  x$four_limb_tissue_valid <- rowSums(x[, tissue_vars, drop = FALSE] != 0 | is.na(x[, tissue_vars, drop = FALSE])) == 0
  x$four_limb_lean_complete <- complete.cases(x[, lean_vars, drop = FALSE])
  x$strict_ALM_valid <- x$DXAEXSTS == 1 & x$four_limb_tissue_valid & x$four_limb_lean_complete
  x$ALM_kg_strict <- ifelse(x$strict_ALM_valid, rowSums(x[, lean_vars, drop = FALSE]) / 1000, NA_real_)
  x$ALMI_strict <- ifelse(!is.na(x$ALM_kg_strict) & !is.na(x$BMXHT) & x$BMXHT > 0,
                          x$ALM_kg_strict / (x$BMXHT / 100)^2, NA_real_)

  x$survey_strata_id <- interaction(x$survey_cycle, x$SDMVSTRA, drop = TRUE, lex.order = TRUE)
  x$survey_psu_id <- interaction(x$survey_cycle, x$SDMVSTRA, x$SDMVPSU, drop = TRUE, lex.order = TRUE)

  x$survey_base_eligible <- x$RIDAGEYR >= 20 & x$RIDAGEYR <= 59 & x$pbcd_available &
    !is.na(x$weight_2yr_raw) & x$weight_2yr_raw > 0 & !is.na(x$SDMVSTRA) & !is.na(x$SDMVPSU)

  x$core_analysis_eligible_strict <- x$survey_base_eligible & !is.na(x$LBXBCD) &
    x$strict_ALM_valid & !is.na(x$BMXHT) & x$BMXHT > 0 & !is.na(x$ALMI_strict)

  x
}

parts <- lapply(seq_len(nrow(cycle_spec)), function(i) build_cycle_core(cycle_spec[i, , drop = FALSE]))
master <- do.call(rbind, parts)
survey_base <- master[master$survey_base_eligible, , drop = FALSE]
strict <- survey_base[survey_base$core_analysis_eligible_strict, , drop = FALSE]

expected_cycle <- c(`2011-2012` = 2604L, `2013-2014` = 1523L, `2015-2016` = 1337L, `2017-2018` = 2164L)
obs_cycle <- table(strict$survey_cycle)
for (cyc in names(expected_cycle)) stop_if_not_equal(obs_cycle[[cyc]], expected_cycle[[cyc]], paste0("Strict N ", cyc))
stop_if_not_equal(nrow(strict), 7628L, "Strict analytic N")
stop_if_not_equal(nrow(survey_base), 10453L, "Survey base N")

saveRDS(survey_base, file.path(generated_dir, "nhanes_survey_base_20_59_core.rds"), compress = "xz")
saveRDS(strict, file.path(generated_dir, "nhanes_core_analysis_20_59_strict.rds"), compress = "xz")

cycle_counts <- data.frame(cycle = names(expected_cycle), strict_core_n = as.integer(obs_cycle[names(expected_cycle)]))
utils::write.csv(cycle_counts, file.path(generated_dir, "core_cycle_counts.csv"), row.names = FALSE)
cat("01_build_core_dataset.R PASS: survey base N=", nrow(survey_base), "; strict N=", nrow(strict), "\n", sep = "")
