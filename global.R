suppressPackageStartupMessages({
  library(shiny)
  library(shinycssloaders)
  library(tidyverse)
  library(ggplot2)
  library(plotly)
  library(pheatmap)
  library(survival)
  library(randomForestSRC)
  library(edgeR)
  library(limma)
  library(DESeq2)
  library(pROC)
})

app_config <- list(
  data_dir = "top2000",
  cache_dir = "cache",
  top_n_genes = 2000,
  use_cache = FALSE
)

cancer_codes <- list(
  "Thyroid" = "thca",
  "Liver" = "lihc",
  "Kidney" = "kirc",
  "Breast" = "brca",
  "Colon" = "coad",
  "Bronchus and Lung" = "luad"
)


clinical_fixed_cols <- c(
  "sample", "patient_id", "tissue_status", "event", "time",
  "age", "gender", "race", "primary_site", "disease_type",
  "sample_type_code"
)

# Expression Matrix
.current_expr_log <- NULL
.current_expr_ct  <- NULL

load_expr_log <- function(ct) {
  if (!identical(ct, .current_expr_ct)) {
    message("Loading expression matrix for: ", ct)
    
    code <- cancer_codes[[ct]]
    file <- file.path(app_config$data_dir,
                      paste0("TCGA_", toupper(code), "_merged_2000genes.csv"))
    
    df <- suppressMessages(readr::read_csv(file, show_col_types = FALSE))
    
    # TRUE gene columns = numeric and not clinical
    gene_cols <- setdiff(
      names(df)[sapply(df, is.numeric)],
      c("age", "time", "sample_type_code")
    )
    
    expr_mat <- as.matrix(t(df[, gene_cols, drop = FALSE]))
    rownames(expr_mat) <- gene_cols
    colnames(expr_mat) <- df$sample
    
    .current_expr_log <<- log2(expr_mat + 1)
    .current_expr_ct  <<- ct
  }
  
  return(.current_expr_log)
}

get_expr_log <- function(ct) load_expr_log(ct)

# Merged Data Loading
.current_merged_df <- NULL
.current_merged_ct <- NULL

load_merged_data <- function(ct) {
  if (!identical(ct, .current_merged_ct)) {
    
    message("Loading merged clinical data for: ", ct)
    code <- cancer_codes[[ct]]
    file <- file.path(app_config$data_dir,
                      paste0("TCGA_", toupper(code), "_merged_2000genes.csv"))
    
    df <- suppressMessages(readr::read_csv(file, show_col_types = FALSE))
    clinical_cols_present <- intersect(clinical_fixed_cols, names(df))
    gene_cols <- setdiff(
      names(df)[sapply(df, is.numeric)],
      c("age", "time", "sample_type_code")
    )
    
    expr_sub <- df[, gene_cols, drop = FALSE]
    vars <- apply(expr_sub, 2, var, na.rm = TRUE)
    
    top_n <- min(app_config$top_n_genes, length(vars))
    top_genes <- names(sort(vars, decreasing = TRUE))[1:top_n]
    
    df_final <- cbind(
      df[, clinical_cols_present, drop = FALSE],
      df[, top_genes, drop = FALSE]
    )
    
    .current_merged_df <<- df_final
    .current_merged_ct <<- ct
  }
  
  return(.current_merged_df)
}


get_merged_data <- function(ct) load_merged_data(ct)


# Survival Data
.current_surv_df <- NULL
.current_surv_ct <- NULL

load_surv_df <- function(ct) {
  if (!identical(ct, .current_surv_ct)) {
    
    message("Loading survival-cleaned data for: ", ct)
    df <- load_merged_data(ct)
    
    df2 <- df %>%
      filter(tissue_status == "Primary Tumor") %>%
      mutate(
        event = case_when(event == "Dead" ~ 1,
                          event == "Alive" ~ 0,
                          TRUE ~ NA_real_),
        time = as.numeric(time),
        age = as.numeric(age),
        gender = as.factor(gender)
      ) %>%
      drop_na(time, event, age, gender)
    
    .current_surv_df <<- df2
    .current_surv_ct <<- ct
  }
  
  return(.current_surv_df)
}

get_surv_df <- function(ct) load_surv_df(ct)

# Gene List
get_gene_list <- function(ct) {
  expr <- get_expr_log(ct)
  
  gene_variance <- apply(expr, 1, var, na.rm = TRUE)
  ranked <- names(sort(gene_variance, decreasing = TRUE))
  
  head(ranked, 10)
}
