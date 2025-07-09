# ============================================================================
# Title       : Primary GLMM Analysis - POCUS Study
# Author      : Ankush Jamthikar, Partho Sengupta
# Description : Load required packages and input dataset for primary analysis
# Date        : Sys.Date()
# ============================================================================

# -----------------------------
# Set working directory
# -----------------------------
setwd("C:/AJ_Contents/(1) Projects/ButterflyPOCUSStudy/(4b) Rscripts/PrimaryAnalysis-GLMM")

# -----------------------------
# Clear workspace
# -----------------------------
rm(list = ls())  # Remove all objects from the environment

# -----------------------------
# Load Required Libraries
# -----------------------------
# Uncomment below lines to install missing packages
# install.packages("readxl")
# install.packages("glmmTMB")


library(readxl)     # To read Excel files
library(irr)        # For calculating Intraclass Correlation Coefficient (ICC)
library(nlme)       # For fitting linear mixed-effects models
library(glmmTMB)    # For fitting generalized linear mixed models
library(lmtest)     # For likelihood ratio tests and model comparison
library(sandwich)   # For robust standard errors (sandwich estimators)
library(survival)   # For survival analysis (if needed later)


# -----------------------------
# Define Input File Path
# -----------------------------
infilepath <- "C:/AJ_Contents/(1) Projects/ButterflyPOCUSStudy/(2) Excelsheet/RevisedData-02-07-2025/BUTTERFLY Mastertracker_v3_208P_06_30_2025.xlsx"

# -----------------------------
# Load the Excel Data
# -----------------------------
in_df <- read_excel(infilepath)

# -----------------------------
# Preview the Loaded Data
# -----------------------------
head(in_df)        # Display first few rows
str(in_df)         # Structure of the dataset
summary(in_df)     # Summary statistics


# ============================================================================
# Step 1: Define Column Mapping (alias → actual column name)
# ============================================================================
col_map <- list(
  los                   = "Length of stay (days)",
  status                = "Discharge",
  intervention          = "Butterfly_1_0",
  intervention_new      = "Butterfly_1_0 withCrossedOvercontrols",
  intervention_sens     = "Butterfly_1_0_sensitivity", 
  clusters              = "team_clusters",
  time                  = "time_clusters",
  time_days             = "time continous (days)",
  time_days_PS          = "time_PS",
  probability           = "GMM_prob0",
  pooled_probability    = "pooled_GMM_prob0",
  cong_card_failure     = "Congestive Cardiac Failure",
  cancer                = "Cancer (type, active/remission, metastasis)",
  Diabetes              = "Diabetes",
  cardiac_arr           = "Cardiac arrhythmia",
  pul_htn               = "Pulmonary hypertension",
  ling_disease          = "lung_diseases_combined",
  Stroke                = "Stroke",
  Venous_thromboembolism = "Venous thromboembolism"
)


# ============================================================================
# Step 2: Filter Dataset for Complete Cases in Key Columns
# ============================================================================
# Select critical columns for initial filtering
cols_to_check <- unlist(col_map[c("los", "intervention", "intervention_new")])

# Keep rows with no missing values in these key columns
df <- in_df[complete.cases(in_df[, cols_to_check]), ]

cat("Filtered Rows:", nrow(df), "\nFiltered Columns:", ncol(df), "\n")


# ============================================================================
# Step 3: Add Clean Aliases as New Columns for Easy Reference
# ============================================================================
# Create new columns in the filtered dataset with simplified alias names
for (alias in names(col_map)) {
  df[[alias]] <- df[[col_map[[alias]]]]
}


# ============================================================================
# Step 4: Compute ICC Using Linear Mixed-Effects Model (LME)
# ============================================================================
# Purpose: Estimate the degree of clustering (ICC) using random intercepts for clusters
# ----------------------------------------------------------------------------

# Define variables for modeling
los          <- df$los                   # Length of stay (response variable)
clusters     <- df$clusters              # Cluster grouping variable
time         <- df$time                  # Categorical time variable (for stepped-wedge)
intervention <- df$intervention          # Binary intervention indicator (0 = control, 1 = POCUS)

# ----------------------------------------------------------------------------
# Fit Linear Mixed-Effects Model
# ----------------------------------------------------------------------------
# Model: log-transformed LOS ~ intervention + time, with random intercepts by cluster
model_lme <- lme(log(los) ~ 1,
                 random = ~1 | clusters)

# Model Summary
summary(model_lme)

# Extract and view variance components
VarCorr(model_lme)

# ----------------------------------------------------------------------------
# Calculate Intraclass Correlation Coefficient (ICC)
# ----------------------------------------------------------------------------
var_components <- VarCorr(model_lme)

# ICC = between-cluster variance / (between-cluster + residual variance)
icc_nlme <- as.numeric(var_components[1,1]) /
  (as.numeric(var_components[1,1]) + as.numeric(var_components[2,1]))

# Display ICC
cat("ICC (nlme):", icc_nlme, "\n")


# ============================================================================
# Step 5: Fit Gamma GLM to Assess Intervention Effect (Component 1)
# ============================================================================
# Purpose: Estimate the multiplicative effect of intervention on LOS
# using Gamma distribution with a log link (suitable for positively skewed data)
# ----------------------------------------------------------------------------

# Fit Gamma GLM
glm_model <- glm(los ~ intervention_new,
                 data   = df,
                 family = Gamma(link = "log"))

# Display model summary
summary(glm_model)


# ============================================================================
# Step 6: Sensitivity Analysis - Gamma GLM Excluding ITT Patients
# ============================================================================
# Purpose: Reassess the effect of intervention on LOS after excluding ITT patients.
# The 'intervention_sens' variable represents a stricter classification of intervention status.
# ----------------------------------------------------------------------------

# Fit Gamma GLM for sensitivity analysis
glm_model <- glm(los ~ intervention_sens,
                 data   = df,
                 family = Gamma(link = "log"))

# Display model summary
summary(glm_model)



# ============================================================================
# Step 7: Association Between Comorbidities and Length of Stay (LOS)
# ============================================================================
# Purpose: Fit Gamma GLMMs to assess how each comorbidity relates to LOS,
# while accounting for clustering (random intercept for team_clusters)
# ----------------------------------------------------------------------------

# Define list of comorbidity variable names (as named in col_map)
comorbidities <- c("Diabetes",
                   "cardiac_arr",               # Cardiac arrhythmia
                   "pul_htn",                   # Pulmonary hypertension
                   "ling_disease",              # Combined lung diseases
                   "Stroke",
                   "Venous_thromboembolism")

# Initialize a results table
results <- data.frame(
  Variable       = character(),
  Estimate       = numeric(),
  `Exp(Estimate)`= numeric(),
  `p-value`      = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each comorbidity and fit Gamma GLMM
for (var in comorbidities) {
  # Dynamically build model formula
  formula <- as.formula(paste0("los ~ ", var, " + (1 | clusters)"))
  
  # Fit Gamma GLMM
  model <- glmmTMB(formula, data = df, family = Gamma(link = "log"))
  
  # Extract coefficient summary
  coef_summary <- summary(model)$coefficients$cond
  
  # Append results
  results <- rbind(results, data.frame(
    Variable        = var,
    Estimate        = coef_summary[2, "Estimate"],
    `Exp(Estimate)` = exp(coef_summary[2, "Estimate"]),
    `p-value`       = coef_summary[2, "Pr(>|z|)"]
  ))
}

# Display results
print(results)


# ============================================================================
# Step 8: Association Between Combined VTE and Lung Disease with LOS
# ============================================================================
# Purpose: Combine Venous Thromboembolism (VTE) and Lung Disease into a single
# binary predictor and assess its association with LOS using a Gamma GLM.
# ----------------------------------------------------------------------------

# Create a new binary variable indicating presence of either VTE or Lung Disease
df$VTE_or_LungDisease <- (df$Venous_thromboembolism == 1 | df$ling_disease == 1)

# Fit Gamma GLM with log link
model_glm <- glm(
  los ~ VTE_or_LungDisease,
  data   = df,
  family = Gamma(link = "log")
)

# Display model summary
summary(model_glm)



















