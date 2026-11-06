# agroforestry-TOF-estimation
R code for estimating Trees Outside Forest (TOF) area across MLRAs in Land Resource Region F using systematic sampling and variance estimation.
# Trees Outside Forests (TOF) Area Estimation — LRR F

R code developed during an agroforestry internship to estimate the area 
of Trees Outside Forests (TOF) across Major Land Resource Areas (MLRAs) 
within Land Resource Region F.

## Overview

This project processes systematic sample data to produce sample-weighted 
TOF area estimates and TOF percentage of eligible area, aggregated per 
MLRA per year and rolled up to LRR F per year.

## Inputs

- Excel files converted to CSV, containing MLRA area, MLRA masked area 
  (forested land with trees), and related fields
- Monte Carlo replicate files (e.g. `montecarlo_wide_lrr_F_2020_mlra_54.csv`) 
  providing per-AOI replicate draws of TOF area for misclassification 
  variance estimation

## Methodology

- Count-based sample weights (N/n) looped across AOI, year, and MLRA
- Variance estimation combines:
  - V1: misclassification variance from Monte Carlo replicates
  - V2: design variance from systematic sampling
- One common forest-mask value per AOI applied across years (team decision 
  for consistency)

## Outputs

Per-MLRA-per-year and LRR-F-per-year:
- Sample-weighted TOF area
- TOF as % of eligible area

## Requirements

R (version X.X+) with the following packages:
- dplyr
- readr
- (list the others you actually use)

## Usage

[brief instructions on how to run — e.g., "Place input CSVs in /data, 
then run `source('tof_estimation.R')`"]
