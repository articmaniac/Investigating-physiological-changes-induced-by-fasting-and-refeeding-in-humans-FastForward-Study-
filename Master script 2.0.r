## =========================
## PACKAGES
## =========================
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(broom.mixed)

## =========================
## HELPERS
## =========================
first_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) NA else y[1]
}

## =========================
## PARAMETERS (KEEP ONLY ONCE)
## =========================
lab_vars <- c(
  "Red Blood Cells (Mio/µl)",
  "Basophils abs. (Tsd./µl)",
  "Eosinophils abs. (Tsd./µl)",
  "Immature Granulocytes abs. (Tsd./µl)",
  "White Blood Cells (Tsd./µl)",
  "Lymphocytes abs. (Tsd./µl)",
  "Monocytes abs. (Tsd./µl)",
  "Neutrophils abs. (Tsd./µl)",
  "Platelets (Tsd./µl)",
  "ESR (mm/h)",
  "CRP (mg/l)",
  "RDW (%)",
  "Hemoglobin (g/dl)",
  "Hematocrit (%)",
  "Alkaline Phosphatase (U/l)",
  "Ferritin (µg/l)",
  "Serum Glucose (mg/dl)",
  "Insulin (mIU/l)",
  "GPT/ALT (U/l)",
  "GOT/AST (U/l)",
  "Gamma-GT (U/l)",
  "GFR (ml/min)",
  "Uric Acid (mg/dl)",
  "PDW (fl)",
  "P-LCR (%)",
  "MCH (pg)",
  "MCHC (g/dl)",
  "MCV (fl)",
  "Monocytes (%)",
  "Neutrophils (%)",
  "Basophils (%)",
  "Eosinophils (%)",
  "Lymphocytes (%)",
  "Immature Granulocytes (%)",
  "Total Cholesterol (mg/dl)",
  "Triglycerides (mg/dl)",
  "HDL-C (mg/dl)",
  "Non-HDL Cholesterol (mg/dl)",
  "LDL-C (mg/dl)",
  "LDL/HDL Ratio (ratio)",
  "TG/HDL Ratio"   # keep in lab_vars if you want it available (even if ANCOVA excludes it)
)

percent_vars <- lab_vars[grepl("%", lab_vars)]

delta_vars <- c(
  "Red Blood Cells (Mio/µl)",
  "Basophils abs. (Tsd./µl)",
  "Eosinophils abs. (Tsd./µl)",
  "Immature Granulocytes abs. (Tsd./µl)",
  "White Blood Cells (Tsd./µl)",
  "Lymphocytes abs. (Tsd./µl)",
  "Monocytes abs. (Tsd./µl)",
  "Neutrophils abs. (Tsd./µl)",
  "Platelets (Tsd./µl)",
  "ESR (mm/h)",
  "CRP (mg/l)",
  "Ferritin (µg/l)",
  "Serum Glucose (mg/dl)",
  "Insulin (mIU/l)",
  "GPT/ALT (U/l)",
  "GOT/AST (U/l)",
  "Gamma-GT (U/l)",
  "GFR (ml/min)",
  "Uric Acid (mg/dl)",
  "Alkaline Phosphatase (U/l)",
  "Total Cholesterol (mg/dl)",
  "Triglycerides (mg/dl)",
  "HDL-C (mg/dl)",
  "Non-HDL Cholesterol (mg/dl)",
  "LDL-C (mg/dl)",
  "TG/HDL Ratio",
  "Hemoglobin (g/dl)",
  "MCH (pg)",
  "MCHC (g/dl)",
  "MCV (fl)",
  "PDW (fl)"
)

## =========================
## READ RAW DATA
## =========================
FF_laboratory_raw <- read_csv2(
  "FF_laboratory_ETL2.csv",
  locale = locale(encoding = "ISO-8859-1")
)

FF_fasting <- read_csv("FF_fasting_duration_and_FR_duration.csv") %>%
  mutate(study_id = as.character(study_id))

FF_metadata <- read_excel("C:/Users/megbi/Desktop/BA/Cohort descriptions/FF_metadata_1501026.xlsx") %>%
  select(-sex) %>%
  mutate(study_id = as.character(study_id))

FF_diet <- read_csv("FF_dietary_data_28-11-2025.csv") %>%
  mutate(study_id = as.character(study_id))

## =========================
## CLEAN / STANDARDIZE CORE LABS (LONG)
## =========================
FF_laboratory <- FF_laboratory_raw %>%
  select(-lab_id, -case_id, -range_values, -birthdate) %>%
  mutate(
    datetime = as.Date(datetime, format = "%d.%m.%Y"),
    study_id = as.character(study_id),
    variable = ifelse(variable == "High-Sensitivity CRP (mg/l)", "CRP (mg/l)", variable),
    value = case_when(
      variable == "CRP (mg/l)" & str_detect(value, "<0\\.5") ~ as.character(0.5 / sqrt(2)),
      variable == "CRP (mg/l)" & str_detect(value, "<0\\.16") ~ as.character(0.16 / sqrt(2)),
      variable == "CRP (mg/l)" & str_detect(value, "^<") ~ {
        num_val <- as.numeric(str_extract(value, "[0-9.]+"))
        as.character(num_val / sqrt(2))
      },
      TRUE ~ value
    ),
    value = as.numeric(value)
  ) %>%
  filter(
    grepl("^[0-9]+$", study_id), #get rid of repeaters 
    timepoint %in% c("T0","T1","T2"), #define Timepoints 
    variable %in% lab_vars
  ) %>%
  mutate(timepoint = factor(timepoint, levels = c("T0","T1","T2"))) #set order 

## =========================
## CLEAN METADATA (ONE ROW PER SUBJECT)
## =========================
FF_metadata_clean <- FF_metadata %>%
  filter(grepl("^[0-9]+$", study_id)) %>%
  group_by(study_id) %>%
  summarise(
    age = first_non_na(age),
    BMI_baseline = first_non_na(BMI_baseline),
    .groups = "drop"
  )

## =========================
## CLEAN DIET (ONE ROW PER SUBJECT)
## =========================
diet_lookup <- FF_diet %>%
  group_by(study_id) %>%
  summarise(FR_dietary = first_non_na(FR_dietary), .groups = "drop") %>%
  mutate(
    FR_dietary = stringr::str_trim(FR_dietary),
    diet_keto = case_when(
      FR_dietary == "KETO" ~ "Keto",
      FR_dietary %in% c("normal", "Diabete Diet") ~ "Non-Keto",
      TRUE ~ NA_character_
    ),
    diet_keto = factor(diet_keto, levels = c("Non-Keto", "Keto"))
  )

## =========================
## MASTER LONG DATASET 
## =========================
FF_laboratory_master <- FF_laboratory %>%
  left_join(FF_fasting, by = "study_id") %>%
  left_join(FF_metadata_clean, by = "study_id") %>%
  left_join(diet_lookup %>% select(study_id, diet_keto), by = "study_id")

## =========================
## TG/HDL RATIO AND ADD TO MASTEr
## =========================
TG_HDL_long <- FF_laboratory_master %>%
  filter(variable %in% c("Triglycerides (mg/dl)", "HDL-C (mg/dl)")) %>%
  group_by(study_id, sex, timepoint, fasting_duration_measured, FR_duration, variable) %>%
  slice_min(datetime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(study_id, sex, timepoint, datetime, fasting_duration_measured, FR_duration,
         age, BMI_baseline, diet_keto, variable, value) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  filter(!is.na(`Triglycerides (mg/dl)`),
         !is.na(`HDL-C (mg/dl)`),
         `HDL-C (mg/dl)` > 0) %>%
  transmute(
    study_id, sex,
    variable = "TG/HDL Ratio",
    timepoint, datetime,
    value = `Triglycerides (mg/dl)` / `HDL-C (mg/dl)`,
    fasting_duration_measured, FR_duration,
    age, BMI_baseline, diet_keto
  )

FF_laboratory_master <- bind_rows(FF_laboratory_master, TG_HDL_long)

## =========================
## DELTA DATASET (WIDE)
## =========================
delta_df <- FF_laboratory_master %>%
  filter(variable %in% delta_vars, timepoint %in% c("T0","T1","T2")) %>%
  group_by(study_id, variable, timepoint) %>%
  slice_min(datetime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(study_id, sex, variable, fasting_duration_measured, FR_duration, age, BMI_baseline, diet_keto),
    names_from = timepoint,
    values_from = value
  ) %>%
  mutate(
    delta_T1_T0 = if_else(!is.na(T1) & !is.na(T0), T1 - T0, NA_real_),
    delta_T2_T0 = if_else(!is.na(T2) & !is.na(T0), T2 - T0, NA_real_),
    delta_T2_T1 = if_else(!is.na(T2) & !is.na(T1), T2 - T1, NA_real_)
  )

## =========================
## ANCOVA DATASET (LONG, baseline-adjusted; EXCLUDES % and optionally excludes TG/HDL)
## =========================
FF_laboratory_clean <- FF_laboratory %>%
  arrange(study_id, variable, timepoint, datetime) %>%
  group_by(study_id, variable, timepoint) %>%
  summarise(
    value = first(value[!is.na(value)]),
    sex = first(sex),
    .groups = "drop"
  )

baseline_df <- FF_laboratory_clean %>%
  filter(timepoint == "T0") %>%
  transmute(study_id, variable, baseline = value)

FF_ancova <- FF_laboratory_clean %>%
  filter(timepoint %in% c("T1","T2")) %>%
  left_join(baseline_df, by = c("study_id","variable")) %>%
  left_join(FF_metadata_clean, by = "study_id") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("T1","T2")),
    sex = factor(sex, levels = c("F","M"))
  ) %>%
  filter(!variable %in% percent_vars) %>%
  filter(variable != "TG/HDL Ratio")   # <- remove ratios from ANCOVA (add more lines if needed)

## =========================
## MHO / MUO classification at baseline (T0)
## =========================

# 0) baseline obese definition (BMI >= 30)
obese_ids <- FF_metadata_clean %>%
  filter(!is.na(BMI_baseline), BMI_baseline >= 30) %>%
  select(study_id, BMI_baseline, age)

# markers used to define baseline abnormalities
mho_markers <- c(
  "Triglycerides (mg/dl)",
  "HDL-C (mg/dl)",
  "CRP (mg/l)",
  "GPT/ALT (U/l)",
  "Uric Acid (mg/dl)"
  # optional later: "Gamma-GT (U/l)"
)

# 1) Baseline labs (T0) in wide format
baseline_T0_wide <- FF_laboratory_clean %>%
  filter(timepoint == "T0", variable %in% mho_markers) %>%
  select(study_id, sex, variable, value) %>%
  pivot_wider(names_from = variable, values_from = value)

# 2) Attach BMI + age
baseline_T0_wide <- baseline_T0_wide %>%
  left_join(obese_ids, by = "study_id")

# 3) Create abnormality flags + abnormality count
baseline_T0_flags <- baseline_T0_wide %>%
  mutate(
    obese = !is.na(BMI_baseline) & BMI_baseline >= 30,

    # lipids
    abn_TG  = !is.na(`Triglycerides (mg/dl)`) & `Triglycerides (mg/dl)` >= 150,
    abn_HDL = case_when(
      sex == "M" ~ !is.na(`HDL-C (mg/dl)`) & `HDL-C (mg/dl)` < 40,
      sex == "F" ~ !is.na(`HDL-C (mg/dl)`) & `HDL-C (mg/dl)` < 50,
      TRUE ~ NA
    ),

    # inflammation
    abn_CRP = !is.na(`CRP (mg/l)`) & `CRP (mg/l)` >= 3,

    # liver enzyme (ALT)
    abn_ALT = case_when(
      sex == "M" ~ !is.na(`GPT/ALT (U/l)`) & `GPT/ALT (U/l)` > 40,
      sex == "F" ~ !is.na(`GPT/ALT (U/l)`) & `GPT/ALT (U/l)` > 35,
      TRUE ~ NA
    ),

    # uric acid
    abn_URIC = case_when(
      sex == "M" ~ !is.na(`Uric Acid (mg/dl)`) & `Uric Acid (mg/dl)` > 7.0,
      sex == "F" ~ !is.na(`Uric Acid (mg/dl)`) & `Uric Acid (mg/dl)` > 6.0,
      TRUE ~ NA
    )
  ) %>%
  mutate(
    # count abnormalities, treating NA as FALSE
    n_abn = rowSums(across(starts_with("abn_"), ~ if_else(is.na(.x), FALSE, .x)))
  )

# 4) classify obese people only
obese_MHO_MUO <- baseline_T0_flags %>%
  filter(obese) %>%
  mutate(
    MUO_status = case_when(
      n_abn >= 2 ~ "MUO",
      n_abn <= 1 ~ "MHO",
      TRUE ~ NA_character_
    ),
    MUO_status = factor(MUO_status, levels = c("MHO", "MUO"))
  )

# 5) your requested output table
mho_status_df <- obese_MHO_MUO %>%
  select(study_id, MUO_status)

mho_status_df

# how many obese and how split?
obese_MHO_MUO %>% count(MUO_status)

# missingness in baseline markers among obese
obese_MHO_MUO %>%
  summarise(across(all_of(mho_markers), ~ sum(is.na(.x)), .names = "missing_{.col}"))

# distribution of abnormality counts
obese_MHO_MUO %>% count(n_abn)

# plot distribution
obese_MHO_MUO %>%
  ggplot(aes(x = n_abn)) +
  geom_histogram(binwidth = 1) +
  theme_bw()

TG_HDL_long %>%
  left_join(mho_status_df, by = "study_id") %>%
  group_by(MUO_status, timepoint) %>%
  summarise(mean_ratio = mean(value, na.rm = TRUE), n = n(), .groups = "drop")

# join status onto baseline long data (T0 only)
baseline_T0_long_mho <- FF_laboratory_clean %>%
  filter(timepoint == "T0") %>%
  left_join(mho_status_df, by = "study_id") %>%
  filter(!is.na(MUO_status))

# choose markers to display
plot_markers <- c("Triglycerides (mg/dl)", "HDL-C (mg/dl)", "CRP (mg/l)",
                  "GPT/ALT (U/l)", "Uric Acid (mg/dl)")

ggplot(
  baseline_T0_long_mho %>% filter(variable %in% plot_markers),
  aes(x = MUO_status, y = value)
) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = NULL, y = "Baseline (T0) value", title = "Baseline differences: MHO vs MUO") +
  theme_bw()

tghdl_mho <- FF_laboratory_master %>%
  filter(variable == "TG/HDL Ratio", timepoint %in% c("T0","T1","T2")) %>%
  left_join(mho_status_df, by = "study_id") %>%
  filter(!is.na(MUO_status))

summary_tghdl <- tghdl_mho %>%
  group_by(MUO_status, timepoint) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n),
    .groups = "drop"
  )

TG_HDL_plot_df <- FF_laboratory_master %>%
  filter(variable == "TG/HDL Ratio",
         timepoint %in% c("T0", "T1", "T2")) %>%
  left_join(mho_status_df, by = "study_id") %>%
  filter(!is.na(MUO_status))

n_per_group <- TG_HDL_plot_df %>%
  group_by(MUO_status, timepoint) %>%
  summarise(
    n_subjects = n_distinct(study_id),
    .groups = "drop"
  )

n_per_group

n_labels_TGHDL <- n_per_group %>%
  left_join(
    summary_tghdl %>%
      select(MUO_status, timepoint, mean, se),
    by = c("MUO_status", "timepoint")
  ) %>%
  mutate(
    label = paste0("n = ", n_subjects),
    y_pos = case_when(
      MUO_status == "MUO" ~ mean + 1.96 * se + 0.05 * mean,  # above CI
      MUO_status == "MHO" ~ -Inf                             # bottom
    ),
    vjust = case_when(
      MUO_status == "MUO" ~ 0,
      MUO_status == "MHO" ~ -0.6
    )
  )

n_labels_TGHDL <- n_per_group %>%
  left_join(
    summary_tghdl %>%
      select(MUO_status, timepoint, mean, se),
    by = c("MUO_status", "timepoint")
  ) %>%
  mutate(
    label = paste0("n = ", n_subjects),
    y_pos = case_when(
      MUO_status == "MUO" ~ mean + 1.96 * se + 0.05 * mean,  # above CI
      MUO_status == "MHO" ~ -Inf                             # bottom
    ),
    vjust = case_when(
      MUO_status == "MUO" ~ 0,
      MUO_status == "MHO" ~ -0.6
    )
  )

ggplot(
  summary_tghdl,
  aes(x = timepoint, y = mean, group = MUO_status, color = MUO_status)
) +
  geom_hline(
    yintercept = 3,
    linetype = "dotted",
    linewidth = 0.8,
    color = "grey40"
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
    width = 0.1
  ) +
  geom_text(
    data = n_labels_TGHDL,
    aes(x = timepoint, y = y_pos, label = label, color = MUO_status, vjust = vjust),
    inherit.aes = FALSE,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    x = "Timepoint",
    y = "TG/HDL ratio (mean ± 95% CI)",
    title = "TG/HDL ratio over time by phenotype",
    subtitle = "Dotted line indicates TG/HDL = 3 (insulin-resistance–associated threshold)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# alfred criticism
# HDL is Not so correlated to baseline, but TG extremly.
# So you need to investigate finally why the ratio is changing, it comes from TG changes, HDL changes, or really both
