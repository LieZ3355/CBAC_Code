# 1. Configuration ------------------------------------------------

exposure_rds <- NULL
input_file <- "Severe.tsv"
output_dir <- "Atherosclerosis,Covid-19 Severe"
outcome_label <- "COVID-19 Severe"
outcome_id <- "COVID19_Severe"
presso_alpha <- 0.05
presso_seed <- 1234
presso_min_simulations <- 10000L

required_packages <- c("data.table", "TwoSampleMR", "MRPRESSO", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Install the missing packages: ", paste(missing_packages, collapse = ", "))
}

if (!exists("exp_dat", inherits = TRUE)) {
  if (!is.null(exposure_rds) && file.exists(exposure_rds)) {
    exp_dat <- readRDS(exposure_rds)
  } else {
    stop("Prepare exp_dat first, or set exposure_rds to a valid RDS file.")
  }
}

exp_dat <- as.data.frame(exp_dat)
required_exposure_cols <- c(
  "SNP", "beta.exposure", "se.exposure", "effect_allele.exposure",
  "other_allele.exposure", "id.exposure"
)
missing_exposure_cols <- setdiff(required_exposure_cols, names(exp_dat))
if (length(missing_exposure_cols) > 0L) {
  stop("exp_dat is missing: ", paste(missing_exposure_cols, collapse = ", "))
}
if (nrow(exp_dat) == 0L) stop("exp_dat has no variants.")
if (length(unique(exp_dat$id.exposure)) != 1L) {
  stop("This script expects a single exposure phenotype in exp_dat.")
}
if (anyNA(exp_dat$SNP) || anyDuplicated(exp_dat$SNP)) {
  stop("exp_dat must have unique, non-missing SNP identifiers.")
}
if (!file.exists(input_file)) stop("Outcome GWAS file not found: ", input_file)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_path <- function(filename) file.path(output_dir, filename)
save_csv <- function(x, filename) {
  utils::write.csv(as.data.frame(x), output_path(filename), row.names = FALSE, na = "")
}

# 2. Read and match outcome variants ------------------------------

outcome_columns <- c(
  "rsid", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_meta_AF",
  "ALT", "REF", "all_inv_var_meta_p", "#CHR", "POS"
)

gwas <- data.table::fread(input_file, select = outcome_columns, showProgress = TRUE)
missing_outcome_cols <- setdiff(outcome_columns, names(gwas))
if (length(missing_outcome_cols) > 0L) {
  stop("Outcome GWAS is missing: ", paste(missing_outcome_cols, collapse = ", "))
}

gwas_snp <- as.character(gwas[["rsid"]])
exposure_snp <- as.character(exp_dat$SNP)
matched_index <- match(exposure_snp, gwas_snp)

# A duplicated SNP identifier may refer to ambiguous GWAS entries.
matched_gwas_ids <- unique(exposure_snp[!is.na(matched_index)])
duplicated_gwas_ids <- unique(gwas_snp[duplicated(gwas_snp)])
ambiguous_ids <- intersect(matched_gwas_ids, duplicated_gwas_ids)
if (length(ambiguous_ids) > 0L) {
  save_csv(data.frame(SNP = ambiguous_ids), "ambiguous_GWAS_SNPs.csv")
  stop("Duplicate matched SNP IDs were found in the outcome GWAS. ",
       "Resolve the entries listed in ambiguous_GWAS_SNPs.csv first.")
}

unmatched_snp <- unique(exposure_snp[is.na(matched_index)])
save_csv(data.frame(SNP = unmatched_snp), "unmatched_exposure_SNPs.csv")
cat("Exposure variants:", length(exposure_snp), "\n")
cat("Matched outcome variants:", sum(!is.na(matched_index)), "\n")
cat("Unmatched variants:", length(unmatched_snp), "\n")
if (all(is.na(matched_index))) stop("No SNPs matched the outcome GWAS.")

pick <- as.data.frame(gwas[matched_index[!is.na(matched_index)], ])
pick$trait <- outcome_label
pick$id <- outcome_id
rm(gwas)
gc(verbose = FALSE)

# 3. Format outcome summary statistics ---------------------------

out_dat <- TwoSampleMR::format_data(
  pick,
  type = "outcome",
  snp_col = "rsid",
  phenotype_col = "trait",
  beta_col = "all_inv_var_meta_beta",
  se_col = "all_inv_var_meta_sebeta",
  eaf_col = "all_meta_AF",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "all_inv_var_meta_p",
  chr_col = "#CHR",
  pos_col = "POS",
  id_col = "id"
)

# Preserve the GWAS effect direction corresponding to the ALT allele.
# Do not negate beta.outcome without a documented allele recoding.
save_csv(out_dat, "outcome_formatted.csv")

# 4. Harmonize exposure and outcome ------------------------------

harmonised_all <- TwoSampleMR::harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat,
  action = 2
)
save_csv(harmonised_all, "harmonised_all_SNPs.csv")

valid_rows <- !is.na(harmonised_all$mr_keep) & harmonised_all$mr_keep
for (field in c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome")) {
  valid_rows <- valid_rows & is.finite(harmonised_all[[field]])
}
valid_rows <- valid_rows &
  harmonised_all$se.exposure > 0 & harmonised_all$se.outcome > 0

dat <- as.data.frame(harmonised_all[valid_rows, , drop = FALSE])
rownames(dat) <- NULL
if (anyDuplicated(dat$SNP)) stop("Duplicate SNPs remained after harmonization.")
if (nrow(dat) == 0L) stop("No valid SNPs remained after harmonization.")
save_csv(dat, "harmonised_MR_SNPs.csv")
cat("Harmonized variants retained for MR:", nrow(dat), "\n")

# 5. Helpers for analysis and plot export ------------------------

mr_methods <- function(n_snp) {
  if (n_snp == 1L) return("mr_wald_ratio")
  if (n_snp == 2L) return("mr_ivw_mre")
  c("mr_egger_regression", "mr_weighted_median", "mr_ivw_mre", "mr_weighted_mode")
}

save_first_plot <- function(plot_list, filename) {
  if (is.list(plot_list) && length(plot_list) > 0L &&
      inherits(plot_list[[1]], "ggplot")) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot_list[[1]],
      device = "pdf", width = 8, height = 8
    )
  }
}

run_mr_analysis <- function(mr_dat, analysis_name) {
  analysis_dir <- output_path(analysis_name)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
  analysis_path <- function(filename) file.path(analysis_dir, filename)
  n_snp <- nrow(mr_dat)
  methods <- mr_methods(n_snp)

  mr_result <- TwoSampleMR::mr(mr_dat, method_list = methods)
  utils::write.csv(mr_result, analysis_path("MR_results.csv"), row.names = FALSE)
  utils::write.csv(
    TwoSampleMR::generate_odds_ratios(mr_result),
    analysis_path("MR_odds_ratios.csv"), row.names = FALSE
  )
  utils::write.csv(mr_dat, analysis_path("MR_SNPs.csv"), row.names = FALSE)

  if (n_snp >= 2L) {
    het_methods <- if (n_snp >= 3L) {
      c("mr_ivw", "mr_egger_regression")
    } else {
      "mr_ivw"
    }
    tryCatch({
      het <- TwoSampleMR::mr_heterogeneity(mr_dat, method_list = het_methods)
      utils::write.csv(het, analysis_path("heterogeneity.csv"), row.names = FALSE)
    }, error = function(e) warning("Heterogeneity analysis: ", conditionMessage(e)))
  }
  if (n_snp >= 3L) {
    tryCatch({
      pleio <- TwoSampleMR::mr_pleiotropy_test(mr_dat)
      utils::write.csv(pleio, analysis_path("Egger_intercept.csv"), row.names = FALSE)
    }, error = function(e) warning("Egger intercept: ", conditionMessage(e)))
  }

  tryCatch({
    save_first_plot(
      TwoSampleMR::mr_scatter_plot(mr_result, mr_dat),
      analysis_path("scatter_plot.pdf")
    )
  }, error = function(e) warning("Scatter plot: ", conditionMessage(e)))

  if (n_snp >= 2L) {
    tryCatch({
      res_single <- TwoSampleMR::mr_singlesnp(
        mr_dat,
        all_method = if (n_snp >= 3L) {
          c("mr_ivw_mre", "mr_egger_regression")
        } else {
          "mr_ivw"
        }
      )
      utils::write.csv(res_single, analysis_path("single_SNP_results.csv"), row.names = FALSE)
      save_first_plot(
        TwoSampleMR::mr_forest_plot(res_single), analysis_path("forest_plot.pdf")
      )
      if (n_snp >= 3L) {
        save_first_plot(
          TwoSampleMR::mr_funnel_plot(res_single), analysis_path("funnel_plot.pdf")
        )
      }
    }, error = function(e) warning("Single-SNP or forest/funnel plot: ",
                                  conditionMessage(e)))
  }

  if (n_snp >= 3L) {
    tryCatch({
      res_loo <- TwoSampleMR::mr_leaveoneout(mr_dat)
      utils::write.csv(res_loo, analysis_path("leave_one_out_results.csv"), row.names = FALSE)
      save_first_plot(
        TwoSampleMR::mr_leaveoneout_plot(res_loo),
        analysis_path("leave_one_out_plot.pdf")
      )
    }, error = function(e) warning("Leave-one-out analysis: ", conditionMessage(e)))
  }
  invisible(mr_result)
}

# 6. Primary MR analysis: all eligible harmonized variants -------

res_raw <- run_mr_analysis(dat, "01_Primary_all_SNPs")

# 7. MR-PRESSO: use the same harmonized data as primary MR -------

presso_data <- as.data.frame(dat[, c(
  "SNP", "beta.outcome", "beta.exposure", "se.outcome", "se.exposure"
), drop = FALSE])
rownames(presso_data) <- NULL

# The original MR-PRESSO implementation requires at least four SNPs.
# More simulations improve the precision of Bonferroni-adjusted outlier tests.
n_snp <- nrow(presso_data)
n_sim <- max(presso_min_simulations, ceiling(n_snp / presso_alpha) + 1L)
presso_result <- NULL
presso_status <- "Not run: fewer than four eligible SNPs."

if (n_snp >= 4L) {
  presso_result <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = presso_data,
      NbDistribution = n_sim,
      SignifThreshold = presso_alpha,
      seed = presso_seed
    ),
    error = function(e) {
      warning("MR-PRESSO failed: ", conditionMessage(e))
      NULL
    }
  )
  presso_status <- if (is.null(presso_result)) "Failed; see warnings." else "Completed"
}

# Store the complete MR-PRESSO object, including its original model estimates.
saveRDS(presso_result, output_path("MR_PRESSO_full_result.rds"))
if (!is.null(presso_result)) {
  writeLines(capture.output(print(presso_result)), output_path("MR_PRESSO_report.txt"))
  save_csv(presso_result$`Main MR results`, "MR_PRESSO_main_MR_results.csv")
}

# MR-PRESSO reports per-SNP Bonferroni-adjusted P values. Values
# such as '<0.001' are treated as upper bounds, conservatively.
parse_presso_p <- function(value) {
  raw <- trimws(as.character(value))
  upper_bound <- grepl("^<", raw)
  numeric_value <- suppressWarnings(as.numeric(sub("^<[[:space:]]*", "", raw)))
  list(value = numeric_value, upper_bound = upper_bound)
}

outlier_table <- data.frame(
  SNP = presso_data$SNP,
  PRESSO_Pvalue = NA_character_,
  PRESSO_significant_outlier = FALSE,
  stringsAsFactors = FALSE
)

global_p <- NA_character_
distortion_p <- NA_character_
outlier_snps <- character(0)

if (!is.null(presso_result)) {
  presso_parts <- presso_result$`MR-PRESSO results`
  global_p <- as.character(presso_parts$`Global Test`$Pvalue)
  distortion_p <- if (!is.null(presso_parts$`Distortion Test`)) {
    as.character(presso_parts$`Distortion Test`$Pvalue)
  } else {
    NA_character_
  }
  outlier_test <- presso_parts$`Outlier Test`
  if (!is.null(outlier_test) && nrow(outlier_test) == n_snp) {
    outlier_table$PRESSO_Pvalue <- as.character(outlier_test$Pvalue)
    parsed_p <- parse_presso_p(outlier_table$PRESSO_Pvalue)
    outlier_table$PRESSO_significant_outlier <- !is.na(parsed_p$value) &
      parsed_p$value <= presso_alpha
    outlier_snps <- outlier_table$SNP[
      outlier_table$PRESSO_significant_outlier
    ]
  } else if (!is.null(outlier_test)) {
    warning("MR-PRESSO outlier test row count does not match the input SNP count. ",
            "No SNPs will be removed automatically.")
    presso_status <- paste(presso_status, "Outlier row-count mismatch.")
  }
}

save_csv(outlier_table, "MR_PRESSO_SNP_outlier_test.csv")
save_csv(data.frame(SNP = outlier_snps), "MR_PRESSO_removed_SNPs.csv")
save_csv(data.frame(
  status = presso_status,
  n_snp_tested = n_snp,
  n_simulations = if (n_snp >= 4L) n_sim else NA_integer_,
  global_pvalue = global_p,
  distortion_pvalue = distortion_p,
  n_outliers = length(outlier_snps),
  stringsAsFactors = FALSE
), "MR_PRESSO_summary.csv")

# 8. Outlier-removed sensitivity analysis ------------------------

if (length(outlier_snps) > 0L) {
  dat_clean <- dat[!dat$SNP %in% outlier_snps, , drop = FALSE]
  if (nrow(dat_clean) >= 2L) {
    rownames(dat_clean) <- NULL
    res_clean <- run_mr_analysis(dat_clean, "02_Sensitivity_outliers_removed")
    cat("MR-PRESSO outliers removed:", paste(outlier_snps, collapse = ", "), "\n")
    cat("Variants remaining:", nrow(dat_clean), "\n")
  } else {
    warning("Outlier removal left fewer than two SNPs; corrected MR was skipped.")
  }
} else {
  cat("MR-PRESSO found no removable outliers, or was not available.\n")
}

# 9. Save session information for reproducibility ----------------

writeLines(capture.output(sessionInfo()), output_path("session_info.txt"))
cat("Analysis finished. Results saved to:", normalizePath(output_dir), "\n")
