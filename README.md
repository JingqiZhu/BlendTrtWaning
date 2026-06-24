# BlendTrtWaning

This repository contains data and code needed to reproduce the results in the paper: Zhu, J., Hemstock, M., Che, Z., Baio, G., & Birnie, R. (2026). Flexible Survival Extrapolation with Blended Hazards: Accounting for Treatment Effect Waning in Health Technology Assessment. *Medical decision making*. https://doi.org/10.1177/0272989X261452264. 

This work presents the blended hazard method as a flexible way to account for treatment effect waning while incorporating external evidence in survival extrapolation. NICE TA366 is used as a demonstrating case study.

![base_case_surv_compare](figures/base%20case/survplot_comparison_24_60_5_5.png)

|                        Method                        | Pembrolizumab | Ipilimumab | Increment |
|:----------------------------------------------------:|:-------------:|:----------:|:---------:|
|        7-year RMST from blended hazard method        |      3.59     |    2.83    |    0.76   |
|        Area under updated 7-year Kaplan-Meier        |      3.61     |    2.84    |    0.77   |
| 7-year RMST from piecewise method in TA366 base case |      2.98     |    2.57    |    0.41   |

## Repository Structure

```text
├── data/                 # Data
├── code/                 # R code
├── docs/                 # Rmd document 
├── figures/              # Survival/hazard/HR plots for base case and sensitivity analysis scenarios
├── tables/               # 7-year RMST tables
├── BlendTrtWaning.Rproj  # Project organisation container
├── README.md
```

## R Shiny App

An additional illustrative Shiny App for sensitivity analysis can be accessed [here](https://jzhu919.shinyapps.io/shinyapp/).