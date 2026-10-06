# 01_hdrs_covariate_sensitivity.R
#
# Sensitivity analysis controlling for baseline depressive symptom severity.
#
# Purpose:
#   Re-run TR versus TNR ROI-to-ROI functional connectivity models after adding
#   baseline HAMD17/HDRS score as a covariate.
#
# Model:
#   FC ~ Group + baseline_HAMD17 + Age + Sex + Education
#
# Notes:
#   - The script assumes Group is coded as 1 = treatment responders, 2 = treatment
#     non-responders, and 3 = healthy controls.
#   - Only Groups 1 and 2 are used here.
#   - Edit the CONFIGURATION section to match local file paths and column names.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(writexl)
})

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

input_file <- "path/to/deidentified_analysis_dataset.xlsx"
sheet_name <- 1
output_dir <- "outputs"

# Column names in the de-identified spreadsheet
subject_col <- "Subject_ID"      # not required for modeling, retained for clarity
group_col <- "Group"             # expected coding: 1 = TR, 2 = TNR, 3 = HC
age_col <- "Age"
sex_col <- "Gender"              # change to "Sex" if that is the column name
education_col <- "Education_Level"
hdrs_col <- "HAMD17"             # baseline HDRS/HAMD score

# TR-vs-TNR ROI pairs reported in the main ROI-to-ROI group-comparison results.
# These names should match the FC column names in the dataset.
tr_tnr_fc_pairs <- c(
  "L_PAG_L_AMY",
  "R_vlTHAL_L_IFJ",
  "HY_R_S2",
  "HY_R_SPL",
  "HY_SMA",
  "vmPFC_L_dACC",
  "L_IFJ_R_dACC",
  "SMA_R_parHG",
  "L_AMY_R_PT",
  "R_CD_R_HP"
)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

bt <- function(x) paste0("`", x, "`")

run_hdrs_sensitivity <- function(data, fc, group_col, hdrs_col, covars) {
  vars_needed <- c(fc, group_col, hdrs_col, covars)
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    return(tibble(
      fc = fc,
      n = NA_integer_,
      beta_group_TNR_vs_TR = NA_real_,
      t_group = NA_real_,
      p_group = NA_real_,
      beta_hdrs = NA_real_,
      p_hdrs = NA_real_,
      note = paste("Missing columns:", paste(missing_vars, collapse = ", "))
    ))
  }

  tmp <- data %>%
    select(all_of(vars_needed)) %>%
    drop_na()

  if (nrow(tmp) < 20) {
    return(tibble(
      fc = fc,
      n = nrow(tmp),
      beta_group_TNR_vs_TR = NA_real_,
      t_group = NA_real_,
      p_group = NA_real_,
      beta_hdrs = NA_real_,
      p_hdrs = NA_real_,
      note = "Too few complete observations"
    ))
  }

  formula_txt <- paste(
    bt(fc), "~", group_col, "+", bt(hdrs_col), "+",
    paste(bt(covars), collapse = " + ")
  )

  fit <- lm(as.formula(formula_txt), data = tmp)
  coefs <- summary(fit)$coefficients

  group_term <- paste0(group_col, "TNR")
  if (!group_term %in% rownames(coefs)) {
    group_term <- grep("TNR", rownames(coefs), value = TRUE)[1]
  }

  tibble(
    fc = fc,
    n = nrow(tmp),
    beta_group_TNR_vs_TR = coefs[group_term, "Estimate"],
    t_group = coefs[group_term, "t value"],
    p_group = coefs[group_term, "Pr(>|t|)"],
    beta_hdrs = coefs[hdrs_col, "Estimate"],
    p_hdrs = coefs[hdrs_col, "Pr(>|t|)"],
    note = NA_character_
  )
}

# -----------------------------------------------------------------------------
# LOAD DATA AND RUN ANALYSIS
# -----------------------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_excel(input_file, sheet = sheet_name) %>%
  mutate(
    Group = factor(.data[[group_col]], levels = c(1, 2, 3), labels = c("TR", "TNR", "HC"))
  )

covars <- c(age_col, sex_col, education_col)

dat_mdd <- dat %>%
  filter(Group %in% c("TR", "TNR")) %>%
  mutate(Group = droplevels(Group))

results <- map_dfr(
  tr_tnr_fc_pairs,
  ~ run_hdrs_sensitivity(dat_mdd, .x, "Group", hdrs_col, covars)
) %>%
  mutate(pFDR_group = p.adjust(p_group, method = "BH")) %>%
  arrange(p_group)

print(results, n = Inf)

write.csv(results, file.path(output_dir, "hdrs_covariate_sensitivity_results.csv"), row.names = FALSE)
writexl::write_xlsx(results, file.path(output_dir, "hdrs_covariate_sensitivity_results.xlsx"))

cat("\nNominally significant after controlling for baseline HDRS:\n")
print(results %>% filter(!is.na(p_group), p_group < 0.05), n = Inf)

cat("\nFDR-significant after controlling for baseline HDRS:\n")
print(results %>% filter(!is.na(pFDR_group), pFDR_group < 0.05), n = Inf)
