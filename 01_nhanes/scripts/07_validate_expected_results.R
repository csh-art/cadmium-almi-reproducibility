source(file.path("01_nhanes","scripts","00_config.R"),encoding="UTF-8")
checks <- list()
readg <- function(name) utils::read.csv(file.path(generated_dir,name),check.names=FALSE)
ms <- readg("model_sample_size_validation.csv")
checks[["Model 3 N"]] <- ms$complete_case_n[ms$model=="M3_Primary"] == 6575
checks[["Model 4 N"]] <- ms$complete_case_n[ms$model=="M4_Clinical"] == 6536
q <- readg("weighted_cadmium_quartile_cutpoints.csv")
checks[["Quartile cutpoints"]] <- all(abs(q$weighted_cutpoint_ug_L-c(.16,.25,.47))<1e-12)
reg <- readg("continuous_cadmium_results_all_models.csv")
checks[["Model 3 beta"]] <- abs(reg$estimate[reg$model=="M3_Primary"] - (-0.159643851702317)) < 5e-8
out <- data.frame(check=names(checks),pass=unlist(checks),row.names=NULL)
utils::write.csv(out,file.path(generated_dir,"validation_summary.csv"),row.names=FALSE)
if(!all(out$pass)) stop("One or more validation anchors failed.",call.=FALSE)
cat("All implemented validation anchors PASS.\n")
