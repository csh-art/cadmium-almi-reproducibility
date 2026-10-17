# ============================================================
# Supplementary Table S5
# Gene-level Mendelian randomization: instruments, combined,
# sex-stratified estimates, and Steiger directionality
#
# Creates:
#   1) Word (.docx)
#   2) PDF (.pdf)
#   3) TXT copy of the exact R code used
#
# Output directory:
# E:/1Data/Cd_Sarcopenia_NHANES/reconstruction_2026_10_14/Supplementary Table
#
# Run in RStudio:
# source(
#   "E:/谷歌浏览器下载/05_generate_Supplementary_Table_S5_MR_results.R",
#   encoding = "UTF-8"
# )
# ============================================================

rm(list = ls())
graphics.off()
invisible(gc())

options(
  stringsAsFactors = FALSE,
  width = 300,
  scipen = 999,
  warn = 1
)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
project_root <- "E:/1Data/Cd_Sarcopenia_NHANES"

m11_root <- file.path(
  project_root,
  "new_workspace",
  "03_results",
  "M11_MR"
)

output_root <- file.path(
  m11_root,
  "03_outputs_figures_tables"
)

output_dir <- file.path(
  project_root,
  "reconstruction_2026_10_14",
  "Supplementary Table"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

script_default <- file.path(
  "E:/谷歌浏览器下载",
  "05_generate_Supplementary_Table_S5_MR_results.R"
)

# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------
user_library <- "D:/RStudio/Rpackages"

if (dir.exists(user_library)) {
  .libPaths(
    unique(
      c(
        user_library,
        .libPaths()
      )
    )
  )
}

required_packages <- c(
  "officer",
  "flextable"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "缺少以下 R 包：\n",
      paste(missing_packages, collapse = "\n"),
      "\n\n请先安装后再运行。本脚本不会自动安装包，也不需要 Rtools。"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------
find_named_file <- function(
  root,
  candidate_names
) {

  all_files <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE
  )

  hits <- all_files[
    basename(all_files) %in% candidate_names &
      !grepl(
        "(^|[/\\\\])backup([/\\\\]|$)",
        all_files,
        ignore.case = TRUE,
        perl = TRUE
      )
  ]

  if (length(hits) == 0L) {
    stop(
      paste0(
        "找不到以下结果文件中的任何一个：\n",
        paste(candidate_names, collapse = "\n"),
        "\n\n搜索目录：\n",
        root
      ),
      call. = FALSE
    )
  }

  # Prefer the first candidate filename in the supplied order.
  for (nm in candidate_names) {
    these <- hits[basename(hits) == nm]

    if (length(these) > 0L) {

      if (length(these) == 1L) {
        return(
          normalizePath(
            these,
            winslash = "/",
            mustWork = TRUE
          )
        )
      }

      md5 <- unname(
        as.character(
          tools::md5sum(these)
        )
      )

      if (length(unique(md5)) > 1L) {
        stop(
          paste0(
            "发现多个内容不同的同名文件：",
            nm,
            "\n",
            paste(
              these,
              md5,
              sep = " | MD5=",
              collapse = "\n"
            )
          ),
          call. = FALSE
        )
      }

      return(
        normalizePath(
          these[
            order(
              nchar(these),
              these
            )
          ][1],
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  stop(
    "未能确定结果文件。",
    call. = FALSE
  )
}

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    sep = "\t"
  )
}

pick_col <- function(
  df,
  candidates,
  required = TRUE,
  label = "column"
) {

  nms <- names(df)
  idx <- match(
    tolower(candidates),
    tolower(nms)
  )

  idx <- idx[
    !is.na(idx)
  ]

  if (length(idx) > 0L) {
    return(
      nms[
        idx[1]
      ]
    )
  }

  if (required) {
    stop(
      paste0(
        "无法识别 ",
        label,
        "。\n\n当前字段：\n",
        paste(
          nms,
          collapse = " | "
        )
      ),
      call. = FALSE
    )
  }

  NA_character_
}

num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

wald_p <- function(beta, se) {

  z <- abs(
    beta / se
  )

  2 * stats::pnorm(
    z,
    lower.tail = FALSE
  )
}

fmt_p <- function(x) {

  ifelse(
    is.na(x) |
      !is.finite(x),
    "NA",
    ifelse(
      x < 0.001,
      formatC(
        x,
        format = "e",
        digits = 2
      ),
      formatC(
        x,
        format = "f",
        digits = 4
      )
    )
  )
}

fmt_num <- function(
  x,
  digits = 4,
  signed = TRUE
) {

  ifelse(
    is.na(x) |
      !is.finite(x),
    "NA",
    if (signed) {
      sprintf(
        paste0(
          "%+.",
          digits,
          "f"
        ),
        x
      )
    } else {
      sprintf(
        paste0(
          "%.",
          digits,
          "f"
        ),
        x
      )
    }
  )
}

detect_current_script <- function() {

  frames <- sys.frames()

  ofiles <- vapply(
    frames,
    function(x) {
      if (!is.null(x$ofile)) {
        as.character(x$ofile)[1]
      } else {
        NA_character_
      }
    },
    character(1)
  )

  ofiles <- ofiles[
    !is.na(ofiles) &
      nzchar(ofiles)
  ]

  if (length(ofiles) > 0L) {
    candidate <- tail(
      ofiles,
      1
    )

    if (file.exists(candidate)) {
      return(
        normalizePath(
          candidate,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  if (file.exists(script_default)) {
    return(
      normalizePath(
        script_default,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  NA_character_
}

quote_ps <- function(x) {

  paste0(
    "'",
    gsub(
      "'",
      "''",
      normalizePath(
        x,
        winslash = "\\",
        mustWork = FALSE
      ),
      fixed = TRUE
    ),
    "'"
  )
}

convert_docx_to_pdf_word <- function(
  docx_file,
  pdf_file
) {

  powershell_candidates <- c(
    Sys.which("powershell"),
    "C:/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"
  )

  powershell_candidates <- powershell_candidates[
    nzchar(powershell_candidates) &
      file.exists(powershell_candidates)
  ]

  if (length(powershell_candidates) == 0L) {
    stop(
      "未检测到 Windows PowerShell。",
      call. = FALSE
    )
  }

  powershell_exe <- normalizePath(
    powershell_candidates[1],
    winslash = "/",
    mustWork = TRUE
  )

  if (file.exists(pdf_file)) {
    unlink(
      pdf_file,
      force = TRUE
    )
  }

  ps1_file <- tempfile(
    pattern = "supp_table_s5_word_export_",
    fileext = ".ps1"
  )

  ps_lines <- c(
    "$ErrorActionPreference = 'Stop'",
    paste0(
      "$docx = ",
      quote_ps(docx_file)
    ),
    paste0(
      "$pdf = ",
      quote_ps(pdf_file)
    ),
    "$word = $null",
    "$doc = $null",
    "try {",
    "  $word = New-Object -ComObject Word.Application",
    "  $word.Visible = $false",
    "  $word.DisplayAlerts = 0",
    "  $doc = $word.Documents.Open($docx, $false, $true)",
    "  $doc.SaveAs([ref]$pdf, [ref]17)",
    "  $doc.Close($false)",
    "  $doc = $null",
    "  $word.Quit()",
    "  $word = $null",
    "  [GC]::Collect()",
    "  [GC]::WaitForPendingFinalizers()",
    "  Write-Output 'WORD_PDF_EXPORT_OK'",
    "}",
    "catch {",
    "  if ($doc -ne $null) { try { $doc.Close($false) } catch {} }",
    "  if ($word -ne $null) { try { $word.Quit() } catch {} }",
    "  Write-Error $_.Exception.Message",
    "  exit 1",
    "}"
  )

  writeLines(
    ps_lines,
    ps1_file,
    useBytes = TRUE
  )

  result <- suppressWarnings(
    system2(
      command = powershell_exe,
      args = c(
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-File",
        shQuote(
          normalizePath(
            ps1_file,
            winslash = "\\",
            mustWork = TRUE
          )
        )
      ),
      stdout = TRUE,
      stderr = TRUE,
      wait = TRUE
    )
  )

  status <- attr(
    result,
    "status"
  )

  if (is.null(status)) {
    status <- 0L
  }

  if (file.exists(ps1_file)) {
    unlink(
      ps1_file,
      force = TRUE
    )
  }

  if (
    status != 0L ||
      !file.exists(pdf_file) ||
      is.na(file.info(pdf_file)$size) ||
      file.info(pdf_file)$size < 1000
  ) {
    stop(
      paste0(
        "Microsoft Word PDF 导出失败。\n\n",
        paste(
          result,
          collapse = "\n"
        )
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ------------------------------------------------------------
# 4. Locate source files
# ------------------------------------------------------------
instrument_file <- find_named_file(
  output_root,
  c(
    "Table_MR_instrument_characteristics_and_results.tsv",
    "Table_MR_final_manuscript_summary.tsv"
  )
)

combined_file <- find_named_file(
  output_root,
  c(
    "Table_MR_Wald_ratio_combined_ALM_final.tsv"
  )
)

sex_file <- find_named_file(
  output_root,
  c(
    "Table_MR_Wald_ratio_sex_stratified_ALM.tsv"
  )
)

steiger_file <- find_named_file(
  output_root,
  c(
    "Table_MR_Steiger_directionality.tsv"
  )
)

instrument_raw <- read_tsv(
  instrument_file
)

combined_raw <- read_tsv(
  combined_file
)

sex_raw <- read_tsv(
  sex_file
)

steiger_raw <- read_tsv(
  steiger_file
)

cat(
  "[INFO] Instrument fields: ",
  paste(
    names(instrument_raw),
    collapse = " | "
  ),
  "\n",
  sep = ""
)

cat(
  "[INFO] Combined MR fields: ",
  paste(
    names(combined_raw),
    collapse = " | "
  ),
  "\n",
  sep = ""
)

cat(
  "[INFO] Sex-stratified MR fields: ",
  paste(
    names(sex_raw),
    collapse = " | "
  ),
  "\n",
  sep = ""
)

cat(
  "[INFO] Steiger fields: ",
  paste(
    names(steiger_raw),
    collapse = " | "
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 5. Locked target genes
# ------------------------------------------------------------
prespecified_genes <- c(
  "ATF4",
  "CASP9",
  "CYCS",
  "CASP3",
  "PPARGC1A",
  "MAP1LC3B"
)

retained_genes <- c(
  "ATF4",
  "CASP9",
  "CASP3",
  "MAP1LC3B"
)

# ------------------------------------------------------------
# 6. Parse instrument table
# ------------------------------------------------------------
ig <- pick_col(
  instrument_raw,
  c(
    "gene_name",
    "Gene",
    "gene"
  ),
  TRUE,
  "instrument gene"
)

irs <- pick_col(
  instrument_raw,
  c(
    "rsid",
    "SNP",
    "snp"
  ),
  TRUE,
  "instrument rsid"
)

ifstat <- pick_col(
  instrument_raw,
  c(
    "F_stat",
    "F",
    "Fstat",
    "F_statistic"
  ),
  TRUE,
  "instrument F statistic"
)

instrument <- data.frame(
  Gene = as.character(
    instrument_raw[[ig]]
  ),
  SNP = as.character(
    instrument_raw[[irs]]
  ),
  F_stat = num(
    instrument_raw[[ifstat]]
  ),
  stringsAsFactors = FALSE
)

instrument <- instrument[
  instrument$Gene %in%
    retained_genes,
  ,
  drop = FALSE
]

instrument <- instrument[
  !duplicated(
    instrument$Gene
  ),
  ,
  drop = FALSE
]

if (!setequal(
  instrument$Gene,
  retained_genes
)) {
  stop(
    paste0(
      "最终工具变量未能恢复出4个预期基因。\n",
      "Observed: ",
      paste(
        instrument$Gene,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 7. Parse combined MR
# ------------------------------------------------------------
cg <- pick_col(
  combined_raw,
  c(
    "gene_name",
    "Gene",
    "gene"
  ),
  TRUE,
  "combined gene"
)

crs <- pick_col(
  combined_raw,
  c(
    "rsid",
    "SNP",
    "snp"
  ),
  TRUE,
  "combined rsid"
)

cb <- pick_col(
  combined_raw,
  c(
    "beta_MR",
    "beta",
    "estimate"
  ),
  TRUE,
  "combined beta"
)

cs <- pick_col(
  combined_raw,
  c(
    "se_MR",
    "SE",
    "se"
  ),
  TRUE,
  "combined SE"
)

cp <- pick_col(
  combined_raw,
  c(
    "p_MR",
    "P_value",
    "pvalue",
    "P",
    "p"
  ),
  FALSE
)

clo <- pick_col(
  combined_raw,
  c(
    "ci_lower",
    "CI_lower",
    "lower95",
    "Lower95"
  ),
  FALSE
)

chi <- pick_col(
  combined_raw,
  c(
    "ci_upper",
    "CI_upper",
    "upper95",
    "Upper95"
  ),
  FALSE
)

combined <- data.frame(
  Gene = as.character(
    combined_raw[[cg]]
  ),
  SNP = as.character(
    combined_raw[[crs]]
  ),
  beta = num(
    combined_raw[[cb]]
  ),
  SE = num(
    combined_raw[[cs]]
  ),
  P = if (!is.na(cp)) {
    num(
      combined_raw[[cp]]
    )
  } else {
    NA_real_
  },
  CI_low = if (!is.na(clo)) {
    num(
      combined_raw[[clo]]
    )
  } else {
    NA_real_
  },
  CI_high = if (!is.na(chi)) {
    num(
      combined_raw[[chi]]
    )
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

combined <- combined[
  combined$Gene %in%
    retained_genes,
  ,
  drop = FALSE
]

combined <- combined[
  !duplicated(
    combined$Gene
  ),
  ,
  drop = FALSE
]

if (nrow(combined) != 4L) {
  stop(
    paste0(
      "Combined MR 应为4个基因，当前为 ",
      nrow(combined),
      "。"
    ),
    call. = FALSE
  )
}

combined$P[
  !is.finite(
    combined$P
  )
] <- wald_p(
  combined$beta[
    !is.finite(
      combined$P
    )
  ],
  combined$SE[
    !is.finite(
      combined$P
    )
  ]
)

combined$CI_low[
  !is.finite(
    combined$CI_low
  )
] <- combined$beta[
  !is.finite(
    combined$CI_low
  )
] - 1.96 *
  combined$SE[
    !is.finite(
      combined$CI_low
    )
  ]

combined$CI_high[
  !is.finite(
    combined$CI_high
  )
] <- combined$beta[
  !is.finite(
    combined$CI_high
  )
] + 1.96 *
  combined$SE[
    !is.finite(
      combined$CI_high
    )
  ]

combined$P_Bonferroni <- pmin(
  combined$P * 4,
  1
)

combined$P_BH <- p.adjust(
  combined$P,
  method = "BH"
)

# ------------------------------------------------------------
# 8. Parse Steiger table
# ------------------------------------------------------------
sg <- pick_col(
  steiger_raw,
  c(
    "gene_name",
    "Gene",
    "gene"
  ),
  TRUE,
  "Steiger gene"
)

srs <- pick_col(
  steiger_raw,
  c(
    "rsid",
    "SNP",
    "snp"
  ),
  TRUE,
  "Steiger rsid"
)

sr2e <- pick_col(
  steiger_raw,
  c(
    "r2_exposure",
    "R2_exposure",
    "r2_expression"
  ),
  TRUE,
  "Steiger expression R2"
)

sr2o <- pick_col(
  steiger_raw,
  c(
    "r2_outcome",
    "R2_outcome",
    "r2_ALM"
  ),
  TRUE,
  "Steiger outcome R2"
)

steiger <- data.frame(
  Gene = as.character(
    steiger_raw[[sg]]
  ),
  SNP = as.character(
    steiger_raw[[srs]]
  ),
  R2_expression = num(
    steiger_raw[[sr2e]]
  ),
  R2_ALM = num(
    steiger_raw[[sr2o]]
  ),
  stringsAsFactors = FALSE
)

steiger <- steiger[
  steiger$Gene %in%
    retained_genes,
  ,
  drop = FALSE
]

steiger <- steiger[
  !duplicated(
    steiger$Gene
  ),
  ,
  drop = FALSE
]

if (nrow(steiger) != 4L) {
  stop(
    "Steiger结果应为4个基因。",
    call. = FALSE
  )
}

steiger$Direction_supported <- ifelse(
  steiger$R2_expression >
    steiger$R2_ALM,
  "Yes",
  "No"
)

# ------------------------------------------------------------
# 9. Panel A: all six prespecified genes
# ------------------------------------------------------------
panel_a_raw <- data.frame(
  Gene = prespecified_genes,
  stringsAsFactors = FALSE
)

panel_a_raw$Instrument_status <- ifelse(
  panel_a_raw$Gene %in%
    retained_genes,
  "Eligible cis-eQTL retained",
  "No eligible instrument"
)

panel_a_raw$SNP <- instrument$SNP[
  match(
    panel_a_raw$Gene,
    instrument$Gene
  )
]

panel_a_raw$F_stat <- instrument$F_stat[
  match(
    panel_a_raw$Gene,
    instrument$Gene
  )
]

panel_a_raw$beta <- combined$beta[
  match(
    panel_a_raw$Gene,
    combined$Gene
  )
]

panel_a_raw$CI_low <- combined$CI_low[
  match(
    panel_a_raw$Gene,
    combined$Gene
  )
]

panel_a_raw$CI_high <- combined$CI_high[
  match(
    panel_a_raw$Gene,
    combined$Gene
  )
]

panel_a_raw$P <- combined$P[
  match(
    panel_a_raw$Gene,
    combined$Gene
  )
]

panel_a_raw$P_Bonferroni <- combined$P_Bonferroni[
  match(
    panel_a_raw$Gene,
    combined$Gene
  )
]

panel_a_raw$P_BH <- combined$P_BH[
  match(
    panel_a_raw$Gene,
    combined$Gene
  )
]

panel_a_raw$R2_expression <- steiger$R2_expression[
  match(
    panel_a_raw$Gene,
    steiger$Gene
  )
]

panel_a_raw$R2_ALM <- steiger$R2_ALM[
  match(
    panel_a_raw$Gene,
    steiger$Gene
  )
]

panel_a_raw$Steiger <- steiger$Direction_supported[
  match(
    panel_a_raw$Gene,
    steiger$Gene
  )
]

panel_a <- data.frame(
  Gene = panel_a_raw$Gene,
  `Instrument status` = panel_a_raw$Instrument_status,
  SNP = ifelse(
    is.na(panel_a_raw$SNP),
    "NA",
    panel_a_raw$SNP
  ),
  `F statistic` = ifelse(
    is.finite(panel_a_raw$F_stat),
    sprintf(
      "%.2f",
      panel_a_raw$F_stat
    ),
    "NA"
  ),
  `Wald ratio beta` = fmt_num(
    panel_a_raw$beta,
    digits = 4
  ),
  `95% CI` = ifelse(
    is.finite(panel_a_raw$CI_low) &
      is.finite(panel_a_raw$CI_high),
    paste0(
      fmt_num(
        panel_a_raw$CI_low,
        digits = 4
      ),
      " to ",
      fmt_num(
        panel_a_raw$CI_high,
        digits = 4
      )
    ),
    "NA"
  ),
  `Nominal P` = fmt_p(
    panel_a_raw$P
  ),
  `Bonferroni P` = fmt_p(
    panel_a_raw$P_Bonferroni
  ),
  `BH-adjusted P` = fmt_p(
    panel_a_raw$P_BH
  ),
  `R2 expression` = ifelse(
    is.finite(
      panel_a_raw$R2_expression
    ),
    formatC(
      panel_a_raw$R2_expression,
      format = "e",
      digits = 3
    ),
    "NA"
  ),
  `R2 ALM` = ifelse(
    is.finite(
      panel_a_raw$R2_ALM
    ),
    formatC(
      panel_a_raw$R2_ALM,
      format = "e",
      digits = 3
    ),
    "NA"
  ),
  `Steiger direction` = ifelse(
    is.na(
      panel_a_raw$Steiger
    ),
    "NA",
    panel_a_raw$Steiger
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10. Panel B: sex-stratified results
# ------------------------------------------------------------
xg <- pick_col(
  sex_raw,
  c(
    "gene_name",
    "Gene",
    "gene"
  ),
  TRUE,
  "sex-stratified gene"
)

xrs <- pick_col(
  sex_raw,
  c(
    "rsid",
    "SNP",
    "snp"
  ),
  TRUE,
  "sex-stratified rsid"
)

xsex <- pick_col(
  sex_raw,
  c(
    "sex_group",
    "Sex",
    "sex"
  ),
  TRUE,
  "sex group"
)

xb <- pick_col(
  sex_raw,
  c(
    "beta_MR",
    "beta",
    "estimate"
  ),
  TRUE,
  "sex beta"
)

xs <- pick_col(
  sex_raw,
  c(
    "se_MR",
    "SE",
    "se"
  ),
  TRUE,
  "sex SE"
)

xp <- pick_col(
  sex_raw,
  c(
    "p_MR",
    "P_value",
    "pvalue",
    "P",
    "p"
  ),
  FALSE
)

xlo <- pick_col(
  sex_raw,
  c(
    "ci_lower",
    "CI_lower",
    "lower95",
    "Lower95"
  ),
  FALSE
)

xhi <- pick_col(
  sex_raw,
  c(
    "ci_upper",
    "CI_upper",
    "upper95",
    "Upper95"
  ),
  FALSE
)

sex <- data.frame(
  Gene = as.character(
    sex_raw[[xg]]
  ),
  SNP = as.character(
    sex_raw[[xrs]]
  ),
  Sex = as.character(
    sex_raw[[xsex]]
  ),
  beta = num(
    sex_raw[[xb]]
  ),
  SE = num(
    sex_raw[[xs]]
  ),
  P = if (!is.na(xp)) {
    num(
      sex_raw[[xp]]
    )
  } else {
    NA_real_
  },
  CI_low = if (!is.na(xlo)) {
    num(
      sex_raw[[xlo]]
    )
  } else {
    NA_real_
  },
  CI_high = if (!is.na(xhi)) {
    num(
      sex_raw[[xhi]]
    )
  } else {
    NA_real_
  },
  stringsAsFactors = FALSE
)

sex <- sex[
  sex$Gene %in%
    retained_genes &
    sex$Sex %in%
    c(
      "Female",
      "Male"
    ),
  ,
  drop = FALSE
]

sex$P[
  !is.finite(
    sex$P
  )
] <- wald_p(
  sex$beta[
    !is.finite(
      sex$P
    )
  ],
  sex$SE[
    !is.finite(
      sex$P
    )
  ]
)

sex$CI_low[
  !is.finite(
    sex$CI_low
  )
] <- sex$beta[
  !is.finite(
    sex$CI_low
  )
] - 1.96 *
  sex$SE[
    !is.finite(
      sex$CI_low
    )
  ]

sex$CI_high[
  !is.finite(
    sex$CI_high
  )
] <- sex$beta[
  !is.finite(
    sex$CI_high
  )
] + 1.96 *
  sex$SE[
    !is.finite(
      sex$CI_high
    )
  ]

if (nrow(sex) != 8L) {
  stop(
    paste0(
      "Sex-stratified results should contain 8 rows; observed ",
      nrow(sex),
      "."
    ),
    call. = FALSE
  )
}

# Formal female-male difference, matching the manuscript description.
heterogeneity_records <- lapply(
  retained_genes,
  function(g) {

    f <- sex[
      sex$Gene == g &
        sex$Sex == "Female",
      ,
      drop = FALSE
    ]

    m <- sex[
      sex$Gene == g &
        sex$Sex == "Male",
      ,
      drop = FALSE
    ]

    if (
      nrow(f) != 1L ||
        nrow(m) != 1L
    ) {
      stop(
        "Female/Male rows are not unique for gene: ",
        g,
        call. = FALSE
      )
    }

    z_diff <- (
      f$beta - m$beta
    ) /
      sqrt(
        f$SE^2 +
          m$SE^2
      )

    p_diff <- 2 *
      stats::pnorm(
        abs(z_diff),
        lower.tail = FALSE
      )

    data.frame(
      Gene = g,
      P_sex_difference = p_diff,
      stringsAsFactors = FALSE
    )
  }
)

heterogeneity <- do.call(
  rbind,
  heterogeneity_records
)

sex$P_sex_difference <- heterogeneity$P_sex_difference[
  match(
    sex$Gene,
    heterogeneity$Gene
  )
]

sex$Gene <- factor(
  sex$Gene,
  levels = retained_genes
)

sex$Sex <- factor(
  sex$Sex,
  levels = c(
    "Female",
    "Male"
  )
)

sex <- sex[
  order(
    sex$Gene,
    sex$Sex
  ),
  ,
  drop = FALSE
]

panel_b <- data.frame(
  Gene = as.character(
    sex$Gene
  ),
  Sex = as.character(
    sex$Sex
  ),
  SNP = sex$SNP,
  `Wald ratio beta` = fmt_num(
    sex$beta,
    digits = 4
  ),
  `95% CI` = paste0(
    fmt_num(
      sex$CI_low,
      digits = 4
    ),
    " to ",
    fmt_num(
      sex$CI_high,
      digits = 4
    )
  ),
  `Nominal P` = fmt_p(
    sex$P
  ),
  `P for female-male difference` = fmt_p(
    sex$P_sex_difference
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 11. Locked-result validation
# ------------------------------------------------------------
atf4 <- combined[
  combined$Gene == "ATF4",
  ,
  drop = FALSE
]

if (
  nrow(atf4) != 1L ||
    abs(
      atf4$beta -
        (-0.0165)
    ) > 0.0010
) {
  stop(
    "ATF4 combined MR estimate does not match the locked manuscript result.",
    call. = FALSE
  )
}

if (
  abs(
    min(
      instrument$F_stat,
      na.rm = TRUE
    ) -
      43.46
  ) > 0.20 ||
    abs(
      max(
        instrument$F_stat,
        na.rm = TRUE
      ) -
      311.65
    ) > 0.50
) {
  stop(
    "Instrument F-statistic range does not match the locked manuscript result.",
    call. = FALSE
  )
}

atf4_sex_pdiff <- heterogeneity$P_sex_difference[
  heterogeneity$Gene == "ATF4"
]

if (
  !is.finite(
    atf4_sex_pdiff
  ) ||
    abs(
      atf4_sex_pdiff -
        0.291
    ) > 0.02
) {
  warning(
    paste0(
      "ATF4 sex-difference P recomputed as ",
      signif(
        atf4_sex_pdiff,
        4
      ),
      "; manuscript reports 0.291. Please inspect if the source uses a more exact formula/input."
    ),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 12. Word formatting
# ------------------------------------------------------------
main_title <- paste0(
  "Supplementary Table S5. Gene-level Mendelian randomization analyses ",
  "of skeletal-muscle gene expression and appendicular lean mass"
)

panel_a_title <- paste0(
  "Panel A. Prespecified genes, retained cis-eQTL instruments, ",
  "combined-sample Wald-ratio estimates, and Steiger directionality"
)

panel_b_title <- paste0(
  "Panel B. Sex-stratified Wald-ratio estimates and female-male heterogeneity tests"
)

notes <- paste0(
  "Notes: Six genes were prespecified for gene-level Mendelian randomization. ",
  "One independent skeletal-muscle cis-eQTL was retained for ATF4, CASP9, CASP3, ",
  "and MAP1LC3B; CYCS and PPARGC1A had no eligible instruments. Combined-sample ",
  "Wald-ratio P values were corrected across the four analyzed genes using both ",
  "Bonferroni and Benjamini-Hochberg procedures. Steiger directionality is marked ",
  "\"Yes\" when the retained instrument explained more variance in skeletal-muscle ",
  "gene expression than in appendicular lean mass. Female-male heterogeneity P values ",
  "were calculated from the difference between sex-specific coefficients divided by ",
  "the square root of the sum of their squared standard errors."
)

abbr <- paste0(
  "Abbreviations: ALM, appendicular lean mass; BH, Benjamini-Hochberg; ",
  "CI, confidence interval; cis-eQTL, cis-expression quantitative trait locus; ",
  "MR, Mendelian randomization; R2, proportion of variance explained."
)

make_ft <- function(
  x,
  gene_col = TRUE,
  font_size = 7.5
) {

  ft <- flextable(
    x
  )

  ft <- theme_booktabs(
    ft
  )

  ft <- font(
    ft,
    fontname = "Times New Roman",
    part = "all"
  )

  ft <- fontsize(
    ft,
    size = font_size,
    part = "all"
  )

  ft <- bold(
    ft,
    part = "header"
  )

  ft <- align(
    ft,
    align = "center",
    part = "header"
  )

  ft <- valign(
    ft,
    valign = "center",
    part = "all"
  )

  ft <- padding(
    ft,
    padding = 2.2,
    part = "all"
  )

  if (
    gene_col &&
      "Gene" %in%
        names(x)
  ) {
    ft <- italic(
      ft,
      j = "Gene",
      italic = TRUE,
      part = "body"
    )
  }

  autofit(
    ft
  )
}

ft_a <- make_ft(
  panel_a,
  TRUE,
  7.3
)

ft_b <- make_ft(
  panel_b,
  TRUE,
  8.1
)

doc <- read_docx()

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      main_title,
      fp_text(
        font.family = "Times New Roman",
        font.size = 12,
        bold = TRUE,
        italic = TRUE
      )
    )
  )
)

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      panel_a_title,
      fp_text(
        font.family = "Times New Roman",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc <- body_add_flextable(
  doc,
  value = ft_a
)

doc <- body_add_par(
  doc,
  ""
)

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      panel_b_title,
      fp_text(
        font.family = "Times New Roman",
        font.size = 10,
        bold = TRUE
      )
    )
  )
)

doc <- body_add_flextable(
  doc,
  value = ft_b
)

doc <- body_add_par(
  doc,
  ""
)

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      notes,
      fp_text(
        font.family = "Times New Roman",
        font.size = 8
      )
    )
  )
)

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      abbr,
      fp_text(
        font.family = "Times New Roman",
        font.size = 8
      )
    )
  )
)

doc <- body_set_default_section(
  doc,
  value = prop_section(
    page_size = page_size(
      orient = "landscape"
    ),
    page_margins = page_mar(
      top = 0.42,
      bottom = 0.42,
      left = 0.42,
      right = 0.42
    )
  )
)

# ------------------------------------------------------------
# 13. Save Word, PDF, and exact code
# ------------------------------------------------------------
docx_file <- file.path(
  output_dir,
  "Supplementary_Table_S5_MR_results.docx"
)

pdf_file <- file.path(
  output_dir,
  "Supplementary_Table_S5_MR_results.pdf"
)

code_txt_file <- file.path(
  output_dir,
  "Code_Supplementary_Table_S5_MR_results.txt"
)

if (file.exists(docx_file)) {
  unlink(
    docx_file,
    force = TRUE
  )
}

print(
  doc,
  target = docx_file
)

if (
  !file.exists(docx_file) ||
    is.na(file.info(docx_file)$size) ||
    file.info(docx_file)$size < 1000
) {
  stop(
    "Word文件生成失败。",
    call. = FALSE
  )
}

convert_docx_to_pdf_word(
  docx_file,
  pdf_file
)

current_script <- detect_current_script()

if (is.na(current_script)) {
  stop(
    "Word和PDF已生成，但无法定位当前R脚本，因此无法保存代码TXT。",
    call. = FALSE
  )
}

writeLines(
  readLines(
    current_script,
    warn = FALSE,
    encoding = "UTF-8"
  ),
  code_txt_file,
  useBytes = TRUE
)

# ------------------------------------------------------------
# 14. Final audit
# ------------------------------------------------------------
files <- c(
  Word = docx_file,
  PDF = pdf_file,
  Code_TXT = code_txt_file
)

audit <- data.frame(
  Type = names(files),
  File = unname(files),
  Exists = file.exists(files),
  Size_KB = round(
    file.info(files)$size / 1024,
    2
  ),
  stringsAsFactors = FALSE
)

if (!all(audit$Exists)) {
  stop(
    "最终输出审计失败。",
    call. = FALSE
  )
}

cat(
  "\n",
  "========== SUPPLEMENTARY TABLE S5 COMPLETE ==========\n",
  "Output directory:\n",
  output_dir,
  "\n\n",
  sep = ""
)

print(
  audit,
  row.names = FALSE
)

cat(
  "\nSource files:\n",
  "Instrument: ",
  instrument_file,
  "\n",
  "Combined MR: ",
  combined_file,
  "\n",
  "Sex-stratified MR: ",
  sex_file,
  "\n",
  "Steiger: ",
  steiger_file,
  "\n\n",
  "Retained instruments: ",
  nrow(instrument),
  "/4\n",
  "F-statistic range: ",
  sprintf(
    "%.2f",
    min(
      instrument$F_stat,
      na.rm = TRUE
    )
  ),
  " to ",
  sprintf(
    "%.2f",
    max(
      instrument$F_stat,
      na.rm = TRUE
    )
  ),
  "\n",
  "ATF4 combined beta: ",
  sprintf(
    "%.4f",
    atf4$beta
  ),
  "\n",
  "ATF4 Bonferroni P: ",
  sprintf(
    "%.4f",
    atf4$P_Bonferroni
  ),
  "\n",
  "ATF4 female-male P difference: ",
  sprintf(
    "%.4f",
    atf4_sex_pdiff
  ),
  "\n",
  "=====================================================\n",
  sep = ""
)
