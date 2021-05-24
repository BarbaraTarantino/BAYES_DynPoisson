# Code and data repository for "Bayesian time-varying autoregressive models of COVID-19 epidemics"

This repository hosts code and data necessary to reproduce the results in: Arkaprava R., Giudici P. and Tarantino B., "Bayesian time-varying autoregressive models of COVID-19 epidemics". Application to Italy and US COVID-19 data is provided. 

## R Code

The R code files allows to fit Bayesian dynamic Poisson log-link models including time-constant or time-varying coefficients for count data. The file structure can be summarised as follows: 

-   `coding_npi.R`: Loads the original COVID-19 data and OxCGRT data, transform ordinal NPI covariates in combination of binary policy measures and creates a wide dataframe that includes count data and policy variables. It creates one dataset for each country specified in the initial character vector. 

-   `main.R`: Main script to fit Bayesian time-costant/time-varying autoregressive models of COVID-19 epidemics. At the end, posterior diagnostics is returned together with the plot of true vs predicted values. 

The `functions` folder contains the main wrappers to fit Bayesian models, obtain posterior diagnostics and some useful plots: 

-   `BAYES_log.R`: The functions to fit time-constant AR(p), INGARCH(p,q), ARX(p), INGARCHX(p,q) model for count data.

-   `tvBAYES_log.R` : The functions to fit time-varying AR(p), INGARCH(p,q), ARX(p), INGARCHX(p,q) model for count data.

-   `HMC.R`: The function to fit Hamiltonian Monte Carlo sampling.

-   `posterior_help.R` and `main_functions.R`: Create posterior diagnostics by summarising informations from model object, i.e. posterior coefficients, predictions and error metrics.

-   `plot.R`: Create descriptive plots from posterior diagnostics object.

## Acknowledgements

This research has received funding from the European Union’s Horizon 2020 research and innovation program “PERISCOPE: Pan European Response to the ImpactS of COvid- 19 and future Pandemics and Epidemics”, under the grant agreement No 101016233, H2020-SC1- PHE CORONAVIRUS-2020-2-RTD).
