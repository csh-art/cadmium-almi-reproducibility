options(stringsAsFactors = FALSE, scipen = 999, survey.lonely.psu = "adjust")

required_packages <- c("haven", "survey")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

find_repo_root <- function() {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidates <- unique(c(wd, dirname(wd), dirname(dirname(wd))))
  hit <- candidates[file.exists(file.path(candidates, "README.md")) & dir.exists(file.path(candidates, "01_nhanes"))]
  if (!length(hit)) stop("Run from the repository root or 01_nhanes/scripts.", call. = FALSE)
  hit[[1]]
}

repo_root <- find_repo_root()
nhanes_dir <- file.path(repo_root, "01_nhanes")
generated_dir <- file.path(nhanes_dir, "generated")
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

raw_xpt_dir <- Sys.getenv("NHANES_RAW_XPT_DIR", unset = file.path(repo_root, "data_raw", "01_raw_xpt"))
if (!dir.exists(raw_xpt_dir)) {
  stop("Raw NHANES folder not found. Set NHANES_RAW_XPT_DIR or create data_raw/01_raw_xpt.", call. = FALSE)
}

cycle_spec <- data.frame(
  cycle = c("2011-2012", "2013-2014", "2015-2016", "2017-2018"),
  suffix = c("G", "H", "I", "J"),
  weight_component = c("DEMO", "PBCD", "PBCD", "DEMO"),
  weight_variable = c("WTMEC2YR", "WTSH2YR", "WTSH2YR", "WTMEC2YR"),
  stringsAsFactors = FALSE
)

read_xpt_checked <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path, call. = FALSE)
  haven::read_xpt(path)
}

make_path <- function(cycle, stem, suffix, special = NULL) {
  name <- if (!is.null(special)) special else paste0(stem, "_", suffix, ".XPT")
  file.path(raw_xpt_dir, cycle, name)
}

stop_if_not_equal <- function(observed, expected, label) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(label, ": observed ", observed, ", expected ", expected, call. = FALSE)
  }
  invisible(TRUE)
}
