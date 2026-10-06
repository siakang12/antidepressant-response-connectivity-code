# 02_group_score_interaction.R
#
# Group-by-score interaction analyses for exploratory brain-behavior associations.
#
# Purpose:
#   Test whether FC-clinical or FC-cognitive association slopes differ between
#   treatment responders and non-responders using a unified regression model.
#
# Model:
#   FC ~ Group * Score_z + Age + Sex + Education
#
# Notes:
#   - This script uses only treatment responders and treatment non-responders.
#   - The score variable is standardized before modeling.
#   - Edit CONFIGURATION and reported_tests if local column names differ.

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
group_col <- "Group"             # expected coding: 1 = TR, 2 = TNR, 3 = HC
age_col <- "Age"
sex_col <- "Gender"              # change to "Sex" if that is the column name
education_col <- "Education_Level"

# Cognitive/clinical score column names. Edit if needed.
score_columns <- list(
  DigitSpan = "DigitSpan",
  MatrixReasoning = "MatrixReasoning",
  VerbalFluency = "Verbalfluency",
  TMT_1 = "TMT_A",
  TMT_2 = "TMT_B",
  StroopTask = "Stroop_color",
  IllnessSeverity = "HAMD17"
)

# FC-score pairs examined in the exploratory correlation analyses.
# FC column names should match the de-identified dataset.
reported_tests <- tibble::tribble(
  ~domain, ~score_key, ~fc,
  "cognition", "DigitSpan",        "vmPFC_L_HP",
  "cognition", "MatrixReasoning",  "L_PAG_L_AMY",
  "cognition", "VerbalFluency",    "L_CD_L_PRE",
  "cognition", "TMT_1",            "R_PAG_R_PT",
  "cognition", "TMT_1",            "R_vlTHAL_R_SMG",
  "cognition", "TMT_2",            "R_vlTHAL_L_IFJ",
  "cognition", "TMT_2",            "L_dACC_L_AMY",
  "cognition", "StroopTask",       "vmPFC_L_HP",
  "cognition", "StroopTask",       "L_dACC_L_AMY",
  "cognition", "StroopTask",       "R_PAG_R_PT",
  "cognition", "StroopTask",       "R_dACC_R_SPL",
  "cognition", "MatrixReasoning",  "L_S2_L_parHG",
  "cognition", "MatrixReasoning",  "R_vlTHAL_HY",
  "cognition", "VerbalFluency",    "R_aINS_R_SPL",
  "cognition", "VerbalFluency",    "R_vlTHAL_R_SMG",
  "cognition", "TMT_1",            "R_aINS_R_SPL",
  "cognition", "TMT_1",            "L_AMY_L_HP",
  "cognition", "TMT_2",            "L_AMY_L_HP",
  "cognition", "TMT_2",            "R_aINS_R_SPL",
  "cognition", "StroopTask",       "R_aINS_R_SPL",
  "cognition", "StroopTask",       "L_AMY_L_HP",
  "clinical",  "IllnessSeverity",  "R_CD_R_PRE",
  "clinical",  "IllnessSeverity",  "R_PAG_R_PT"
) %>% distinct(fc, score_key, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

bt <- function(x) paste0("`", x, "`")

run_interaction_model <- function(data, fc, score_col, group_col, covars) {
  vars_needed <- c(fc, score_col, group_col, covars)
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    return(tibble(
      fc = fc,
      score = score_col,
      n = NA_integer_,
      beta_interaction = NA_real_,
      t_interaction = NA_real_,
      p_interaction = NA_real_,
      note = paste("Missing columns:", paste(missing_vars, collapse = ", "))
    ))
  }

  tmp <- data %>%
    select(all_of(vars_needed)) %>%
    drop_na() %>%
    mutate(Score_z = as.numeric(scale(.data[[score_col]])))

  if (nrow(tmp) < 20) {
    return(tibble(
      fc = fc,
      score = score_col,
      n = nrow(tmp),
      beta_interaction = NA_real_,
      t_interaction = NA_real_,
      p_interaction = NA_real_,
      note = "Too few complete observations"
    ))
  }

  formula_txt <- paste(
    bt(fc), "~", group_col, "* Score_z +",
    paste(bt(covars), collapse = " + ")
  )

  fit <- lm(as.formula(formula_txt), data = tmp)
  coefs <- summary(fit)$coefficients
  interaction_terms <- grep("Group.*:Score_z|Score_z:Group", rownames(coefs), value = TRUE)

  if (length(interaction_terms) == 0) {
    return(tibble(
      fc = fc,
      score = score_col,
      n = nrow(tmp),
      beta_interaction = NA_real_,
      t_interaction = NA_real_,
      p_interaction = NA_real_,
      note = "Interaction term not found"
    ))
  }

  term <- interaction_terms[1]
  tibble(
    fc = fc,
    score = score_col,
    n = nrow(tmp),
    beta_interaction = coefs[term, "Estimate"],
    t_interaction = coefs[term, "t value"],
    p_interaction = coefs[term, "Pr(>|t|)"],
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

dat_mdd <- dat %>%
  filter(Group %in% c("TR", "TNR")) %>%
  mutate(Group = droplevels(Group))

covars <- c(age_col, sex_col, education_col)

analysis_pairs <- reported_tests %>%
  mutate(score_col = unname(unlist(score_columns[score_key])))

results <- pmap_dfr(
  list(analysis_pairs$fc, analysis_pairs$score_col),
  ~ run_interaction_model(dat_mdd, ..1, ..2, "Group", covars)
) %>%
  left_join(analysis_pairs, by = c("fc", "score" = "score_col")) %>%
  mutate(pFDR_interaction = p.adjust(p_interaction, method = "BH")) %>%
  arrange(p_interaction)

print(results, n = Inf)

write.csv(results, file.path(output_dir, "group_by_score_interaction_results.csv"), row.names = FALSE)
writexl::write_xlsx(results, file.path(output_dir, "group_by_score_interaction_results.xlsx"))

cat("\nNominally significant interaction effects:\n")
print(results %>% filter(!is.na(p_interaction), p_interaction < 0.05), n = Inf)

cat("\nFDR-significant interaction effects:\n")
print(results %>% filter(!is.na(pFDR_interaction), pFDR_interaction < 0.05), n = Inf)
