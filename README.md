# BasicGSCA_Prime

## Version 1.1.0

### Author:
Gyeongcheol Cho

## Description:
- The **BasicGSCA_Prime** package enables users to estimate and evaluate basic GSCA models.

## Features:
- Estimate GSCA model parameters and calculate their standard errors (SE) along with 95% confidence intervals (CI).
- Assess model performance based on both explanatory and predictive power.
- Compute the PET (Predictor Exclusion Threshold) statistic to evaluate the predictive power of individual predictor components.
- Enable parallel computing for bootstrap sampling.

## Installation:
To use this package in MATLAB:
1. Clone or download the repository:
   ```bash
   git clone https://github.com/GyeongcheolCho/BasicGSCA_Prime.git
   ```
2. Add the package to your MATLAB path:
   ```matlab
    addpath(genpath('BasicGSCA_Prime'))
   ```

## Usage:
- For examples on how to use the package, refer to the `Run_Example_BasicGSCA.m` file. This file demonstrates the implementation of `BasicGSCA()` using the ACSI dataset.

## Compatibility:
- Tested on MATLAB R2023b.
- Likely compatible with earlier MATLAB versions.

### Citation (APA):
- If you use **BasicGSCA_Prime** in your research or publications, please cite it in APA format as follows:

```plaintext
Cho, G. (2024). BasicGSCA_Prime: A package for basic generalized structured component analysis (Version 1.1.0) [Computer software]. GitHub. https://github.com/GyeongcheolCho/BasicGSCA_Prime
```
