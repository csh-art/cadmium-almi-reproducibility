run_ssgsea <- function(
        expression_matrix,
        gene_sets
) {
    
    current_gsva_version <- utils::packageVersion(
        "GSVA"
    )
    
    if (
        current_gsva_version >=
        utils::package_version(
            "1.50.0"
        )
    ) {
        
        parameter_object <- GSVA::ssgseaParam(
            exprData = expression_matrix,
            geneSets = gene_sets,
            normalize = TRUE
        )
        
        score_matrix <- GSVA::gsva(
            parameter_object,
            verbose = TRUE
        )
        
    } else {
        
        score_matrix <- GSVA::gsva(
            expression_matrix,
            gene_sets,
            method = "ssgsea",
            kcdf = "Gaussian",
            ssgsea.norm = TRUE,
            verbose = TRUE
        )
    }
    
    as.matrix(
        score_matrix
    )
}
