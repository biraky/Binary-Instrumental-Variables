# Things to read:

Generalized moments model

LATE estimate (Wald estimator)



# Simulation Strategy

* Instrumental Variable Z: Binary

* Exposure X: Binary

* Confounder U: Binary = correlated with X and Y (doesn't matter)

* Outcome Y: Binary

## Fit: 
- OLS 
- 2sLS
- GMM
- LATE estimator
expect similar results from 2sLS | GMM | LATE

## Estimate:

Bias in covariate effects

RMSE | R^2

## Vary: 
- sample sizes: n = 50, 100, 500, 1000, 5000, 10000
- Cor(Z,X) = Very weak:: rho = 0.25 | Weak:: rho = 0.35 | Moderate:: rho = 0.5 | Strong:: rho = 0.7 | Very strong:: rho = 0.9
