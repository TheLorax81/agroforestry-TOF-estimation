#============================================================================
  # Trees Outside Forest (TOF) Area Estimation for LRR F
  #
  # Goal: estimate the total area of TOF (trees outside forest) within each
  # MLRA in LRR F, for each survey year, plus an LRR-F total per year,
  # using the same systematic-sampling design and two-part variance approach
  # (V1 = misclassification error, V2 = sampling design error) similar to that
  # of LandSampleEstimation_Key.R.
  #
  # Sampling design recap (from LandSampleEstimation_Key.R):
  #   N      = number of population units (grid cells) in the sampling frame
  #   n      = number of units actually drawn into the systematic sample
  #   a      = N / n   (the sampling interval)
  #   sy_wt  = N / n   (expansion weight applied to each sampled unit)
  #   V2     = N^2/n * (1 - 1/a) * var(per-unit sample values)
  #   V1     = var(Random Monte Carlo replicate totals), where each replicate
  #            reclassifies some sampled area between classes to represent
  #            classification/measurement error
  #
  # In the Key script the population units were 9 grid-cell "rows" of a
  # 6-column raster, so each row contributed 6 cells (rowSums). Here, each
  # AOI *is* a single 1 km^2 population unit (an actual systematic sample
  # point), so there is no extra "columns per row" dimension to sum over --
  # each AOI's own value stands in for what rowSums() computed in the
  # original script.
  # ============================================================================

# ---------------------------------------------------------------------------
# 0. PLACEHOLDER CONSTANTS -- replace these when real numbers are available
# ---------------------------------------------------------------------------
# These stand in for the Key script's hardcoded 5%/15% misclassification
# rates. They are used ONLY for MLRA/year combinations that do NOT yet have
# an empirical Monte Carlo file like montecarlo_wide_lrr_F_2020_mlra_54.csv.
# Swap in a real confusion-matrix-derived rate (or better, an empirical MC
# file per MLRA/year) as soon as one exists.
PLACEHOLDER_MISCLASS_RATE_TOF_TO_MASK  <- 0.15  # fraction of TOF area misclassified as mask (forest)
PLACEHOLDER_MISCLASS_RATE_MASK_TO_TOF  <- 0.15  # fraction of mask (forest) area misclassified as TOF
PLACEHOLDER_N_MC_REPS                  <- 1000  # number of Monte Carlo replicates (Key script used 1000)
GRID_CELL_M2                           <- 1000000  # nominal AOI size: 1 km^2, mirrors "each grid cell = 1 ha" in the Key script

#  ---- Common-across-years mask setting (team decision) --------------------
# Use ONE forest-mask value per AOI across all years instead of a
# year-by-year value, so every year's TOF % uses the same denominator.
# MASK_YEAR_MODE options:
#   "average" -- per-AOI mean across the three years (default)
#   2012 / 2016 / 2020 -- pin to that specific year's mask
MASK_YEAR_MODE <- "average"

set.seed(1)  # for reproducibility of the PLACEHOLDER Monte Carlo only, good for
# consistency.

# ---------------------------------------------------------------------------
# 1. LOAD INPUTS
# ---------------------------------------------------------------------------
# All four files are CSVs. aoi_tof_lrr_F.csv is the CSV export of the
# 'aoi_tof' sheet from placeholder_tof_lrr_F.xlsx (one row per AOI per year).
aoi_tof  <- read.csv("placeholder_tof_lrr_F.csv",        stringsAsFactors = FALSE)
mlra_pop <- read.csv("mlraAreas_lrr_F.csv",       stringsAsFactors = FALSE)  # population (census) areas
aoi_smpl <- read.csv("aoiAreasByMlra_lrr_F.csv",  stringsAsFactors = FALSE)  # sample-level summary (cross-check only)

# ---- Freeze the mask across years -----------------------------------------
# Build a per-(aoi_id, mlra_id) frozen mask lookup based on MASK_YEAR_MODE,
# then overwrite each row's mask_area_m2 (and recompute eligible_area_m2)
# so the downstream loop uses the common value regardless of year.
if (MASK_YEAR_MODE == "average") {
  frozen_mask <- aggregate(mask_area_m2 ~ aoi_id + mlra_id,
                           data = aoi_tof, FUN = mean)
} else {
  frozen_mask <- aoi_tof[aoi_tof$year == MASK_YEAR_MODE,
                         c("aoi_id", "mlra_id", "mask_area_m2")]
}
names(frozen_mask)[names(frozen_mask) == "mask_area_m2"] <- "frozen_mask_m2"

aoi_tof <- merge(aoi_tof, frozen_mask, by = c("aoi_id", "mlra_id"), all.x = TRUE)
aoi_tof$mask_area_m2     <- aoi_tof$frozen_mask_m2
aoi_tof$eligible_area_m2 <- aoi_tof$aoi_area_m2 - aoi_tof$mask_area_m2
aoi_tof$frozen_mask_m2   <- NULL

cat("Mask frozen across years using MASK_YEAR_MODE =", MASK_YEAR_MODE, "\n")

# Random Monte Carlo replicates -- currently only for MLRA 60 (symbol
# "54") / 2020. Wide format: one row per aoi_id, one column per replicate.
# NOTE: source data had 50,000 iterations but the file was truncated to
# 16,383 (Excel's column limit). Still plenty for a stable variance --
# if a fuller export arrives later, just drop it in under the same name.
mc_mlra60_2020 <- read.csv("montecarlo_wide_lrr_F_2020_mlra_54.csv", stringsAsFactors = FALSE)

# Registry of which MLRA/year combos have an empirical MC file, mapped to
# the loaded data frame. Add a line here each time a new one arrives.
RANDOM_MC_REGISTRY <- list(
  list(mlra_id = 60, year = 2020, data = mc_mlra60_2020)
)

# ---------------------------------------------------------------------------
# 2. HELPER: population segment count N for a given MLRA
# ---------------------------------------------------------------------------
# N = total number of 1 km^2 cells that could have been sampled in this
# MLRA (the population size). MLRA land area doesn't change year to year,
# so any year's row gives the same N.
get_population_N <- function(mlra_id_val) {
  # MLRA_ID reads in as character because the file also contains an "F"
  # (LRR-level) row, so compare as character explicitly rather than relying
  # on R's automatic coercion.
  pop_row <- mlra_pop[mlra_pop$level == "MLRA" & mlra_pop$MLRA_ID == as.character(mlra_id_val), ][1, ]
  pop_row$total_m2 / GRID_CELL_M2
}

# ---------------------------------------------------------------------------
# 3. HELPER: V1 via Random Monte Carlo replicates (preferred, when available)
# ---------------------------------------------------------------------------
# Random Monte Carlo: computes V1 directly from pre-simulated
# replicate values, one per AOI per rep.
compute_V1_random <- function(mc_wide, sample_rows, sy_wt) {
  rep_cols <- grep("^rep_", names(mc_wide), value = TRUE)
  
  # Align to the AOIs actually in this MLRA/year sample, in case ordering
  # or membership differs between the two files.
  merged <- merge(sample_rows["aoi_id"], mc_wide, by = "aoi_id", all.x = TRUE)
  rep_matrix <- as.matrix(merged[, rep_cols, drop = FALSE])
  storage.mode(rep_matrix) <- "numeric"
  
  # One weighted MLRA-total TOF estimate per replicate...
  mc_totals <- sy_wt * colSums(rep_matrix, na.rm = TRUE)
  # ...and V1 is the variance of that distribution of totals.
  var(mc_totals)
}

# ---------------------------------------------------------------------------
# 4. HELPER: V1 via placeholder fixed-rate Monte Carlo (fallback)
# ---------------------------------------------------------------------------
# Placeholder Monte Carlo: randomly swaps each AOI's TOF and mask areas
# at the placeholder rate to simulate misclassification.
# Replace with real accuracy data when available.
compute_V1_placeholder <- function(sample_rows, sy_wt, n_reps = PLACEHOLDER_N_MC_REPS) {
  n_aoi <- nrow(sample_rows)
  mc_totals <- numeric(n_reps)
  
  for (p in 1:n_reps) {
    tof_reclass  <- sample_rows$tof_area_m2
    mask_reclass <- sample_rows$mask_area_m2
    
    flip_tof_to_mask <- runif(n_aoi) <= PLACEHOLDER_MISCLASS_RATE_TOF_TO_MASK
    flip_mask_to_tof <- runif(n_aoi) <= PLACEHOLDER_MISCLASS_RATE_MASK_TO_TOF
    
    moved_to_mask <- tof_reclass * flip_tof_to_mask
    moved_to_tof  <- mask_reclass * flip_mask_to_tof
    
    tof_reclass <- tof_reclass - moved_to_mask + moved_to_tof
    
    mc_totals[p] <- sy_wt * sum(tof_reclass)
  }
  
  var(mc_totals)
}

# ---------------------------------------------------------------------------
# 5. LOOP OVER AOI, YEAR, AND MLRA to build sample weights + point estimates
# ---------------------------------------------------------------------------
mlra_year_combos <- unique(aoi_tof[, c("mlra_id", "year")])
mlra_year_combos <- mlra_year_combos[order(mlra_year_combos$mlra_id, mlra_year_combos$year), ]

mlra_results <- data.frame()

for (i in seq_len(nrow(mlra_year_combos))) {
  mlra_id_val <- mlra_year_combos$mlra_id[i]
  year_val    <- mlra_year_combos$year[i]
  
  # Get all AOIs sampled in this MLRA for this year.
  # (Boundary AOIs appear as two rows -- one per MLRA -- and each row is
  # treated as its own sample point.)
  sample_rows <- aoi_tof[aoi_tof$mlra_id == mlra_id_val & aoi_tof$year == year_val, ]
  
  N     <- get_population_N(mlra_id_val)
  n     <- nrow(sample_rows)
  a     <- N / n
  sy_wt <- N / n
  
  # ---- Point estimates (weighted totals), mirrors total.forestland.sample ----
  total_tof_est     <- sy_wt * sum(sample_rows$tof_area_m2)
  total_mask_est    <- sy_wt * sum(sample_rows$mask_area_m2)
  total_eligible_est <- sy_wt * sum(sample_rows$eligible_area_m2)
  tof_pct_est <- 100 * total_tof_est / total_eligible_est
  
  # ---- V2: design (systematic sampling) variance component ----
  # Mirrors: V2.forestland <- N^2/n*(1-1/a)*var(rowSums(land[samp,]==20))
  # Here each AOI's tof_area_m2 plays the role rowSums() played per sampled row.
  V2_tof <- N^2 / n * (1 - 1/a) * var(sample_rows$tof_area_m2)
  
  # ---- V1: misclassification variance component ----
  has_random_mc <- FALSE
  for (reg in RANDOM_MC_REGISTRY) {
    if (reg$mlra_id == mlra_id_val && reg$year == year_val) {
      V1_tof     <- compute_V1_random(reg$data, sample_rows, sy_wt)
      V1_method  <- "random_MC"
      has_random_mc <- TRUE
      break
    }
  }
  if (!has_random_mc) {
    V1_tof    <- compute_V1_placeholder(sample_rows, sy_wt)
    V1_method <- "placeholder_MC"
  }
  
  total_var_tof <- V1_tof + V2_tof
  
  mlra_results <- rbind(mlra_results, data.frame(
    level               = "MLRA",
    mlra_id             = mlra_id_val,
    mlra_symbol         = sample_rows$mlra_symbol[1],
    mlra_name           = sample_rows$mlra_name[1],
    year                = year_val,
    N                   = N,
    n                   = n,
    a                   = a,
    sy_wt               = sy_wt,
    total_tof_m2        = total_tof_est,
    total_mask_m2       = total_mask_est,
    total_eligible_m2   = total_eligible_est,
    tof_pct_of_eligible = tof_pct_est,
    V1_tof              = V1_tof,
    V1_method           = V1_method,
    V2_tof              = V2_tof,
    total_var_tof       = total_var_tof,
    se_tof              = sqrt(total_var_tof),
    stringsAsFactors    = FALSE
  ))
}

# ---------------------------------------------------------------------------
# 6. LRR F TOTAL PER YEAR
# ---------------------------------------------------------------------------
# Each MLRA is an independent stratum of the same systematic sample, so the
# LRR-wide total and its variance are just the sums of the per-MLRA totals
# and variances (standard result for independent stratified estimators).
lrr_results <- data.frame()

for (yr in sort(unique(mlra_results$year))) {
  yr_rows <- mlra_results[mlra_results$year == yr, ]
  
  total_tof_lrr      <- sum(yr_rows$total_tof_m2)
  total_eligible_lrr <- sum(yr_rows$total_eligible_m2)
  total_var_lrr       <- sum(yr_rows$total_var_tof)
  
  lrr_results <- rbind(lrr_results, data.frame(
    level               = "LRR",
    mlra_id             = "F",
    mlra_symbol         = "F",
    mlra_name           = "LRR F",
    year                = yr,
    N                   = sum(yr_rows$N),
    n                   = sum(yr_rows$n),
    a                   = NA,
    sy_wt               = NA,
    total_tof_m2        = total_tof_lrr,
    total_mask_m2       = sum(yr_rows$total_mask_m2),
    total_eligible_m2   = total_eligible_lrr,
    tof_pct_of_eligible = 100 * total_tof_lrr / total_eligible_lrr,
    V1_tof              = sum(yr_rows$V1_tof),
    V1_method           = "sum_of_MLRA_components",
    V2_tof              = sum(yr_rows$V2_tof),
    total_var_tof       = total_var_lrr,
    se_tof              = sqrt(total_var_lrr),
    stringsAsFactors    = FALSE
  ))
}

# ---------------------------------------------------------------------------
# 7. COMBINE AND WRITE FINAL RESULT
# ---------------------------------------------------------------------------
final_results <- rbind(mlra_results, lrr_results)
final_results <- final_results[order(final_results$year,
                                     final_results$level,
                                     final_results$mlra_id), ]

write.csv(final_results, "tof_estimates_lrr_F.csv", row.names = FALSE)

cat("Done. Wrote", nrow(final_results), "rows to tof_estimates_lrr_F.csv\n")
print(final_results[, c("level","mlra_symbol","year","total_tof_m2",
                        "tof_pct_of_eligible","V1_method","se_tof")])

