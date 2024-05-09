# Replication package for Co-factor analysis of citation networks

One compelling use of citation networks is to characterize papers by their relationships to the surrounding literature. We propose a method to characterize papers by embedding them into two distinct co-factor spaces: one describing how papers send citations, and the other describing how papers receive citations. This approach presents several challenges. First, older documents cannot cite newer documents, and thus it is not clear that co-factors are even identifiable. We resolve this challenge by developing a co-factor model for asymmetric adjacency matrices with missing lower triangles and showing that identification is possible. We then frame estimation as a matrix completion problem and develop a specialized implementation of matrix completion because prior implementations are memory bound in our setting. Simulations show that our estimator has promising finite sample properties, and that naive approaches fail to recover latent co-factor structure. We leverage our estimator to investigate 255,780 papers published in statistics journals from 1898 to 2024, resulting in the most comprehensive topic model of the statistics literature to date. We find interpretable co-factors corresponding to many statistical subfields, including time series, variable selection, spatial methods, graphical models, GLM(M)s, causal inference, multiple testing, quantile regression, semi-parametrics, dimension reduction, and several more.

## To replicate our computational results

We use [`renv`](https://rstudio.github.io/renv/) to record package dependencies and [`targets`](https://books.ropensci.org/targets/) to coordinate our simulation study and data analysis.

To replicate our results, clone this Github repository to your local computer. Note that this will install a development version of `fastadi` from source, including compiled code, so you will need a functioning C++ toolchain. See installation instructions for `Rcpp` and `RcppArmadillo` if you are running into issues.

Once you have the repository cloned locally, re-create the project library by calling

``` r
# install.packages("renv")
renv::restore()
```

At this point, you should be ready to replicate our simulation and performance comparison results. Our data analysis relies on proprietary data that we cannot release publicly, so our data analysis results cannot be replicated. Nonetheless, we include our data analysis code here for transparency.

The computational work is organized into three distinct `targets` projects (see [here](https://books.ropensci.org/targets/projects.html#multiple-projects)) for details. You will need to build these one at a time.

This is a computationally intensive project. It's likely that, at some point, the targets pipeline will crash for some sundry computational reason. If this is the case, simply re-run `tar_make()`. `tar_make()` will only attempt to re-run incomplete portions of the build pipeline.

``` r
Sys.setenv(TAR_PROJECT = "simulations")
tar_make()

Sys.setenv(TAR_PROJECT = "analysis")
tar_make()

Sys.setenv(TAR_PROJECT = "speed")
tar_make()
```
## Contents

```r
.
├── R
│   ├── analysis-create-graph.R      --  Construct tidygraph object from WoS data extract
│   ├── analysis-estimators.R        --  Estimators for data analysis
│   ├── analysis-forward-citation.R  --  Postprocess imputed forward citations
│   ├── analysis-interpret-vsp.R     --  Create hub tables and plot diagnostics
│   ├── analysis-keywords.R          --  Create keyword tables
│   ├── simulation-estimators.R      --  Estimators for simulation
│   ├── simulation-loss.R            --  Loss calculations
│   ├── simulation-models.R          --  Simulation models
│   └── simulation-vsp.R             --  Align true and estimated factors
├── README.md                        
├── _analysis                        --  [targets]  
├── _analysis.R                      --  Define data analysis pipeline
├── _simulations                     --  [targets]
├── _simulations.R                   --  Define simulation pipeline
├── _speed                           --  [targets]
├── _speed.R                         --  Define performance comparison pipeline
├── _targets.yaml                    --  [targets]
├── citation-cofactoring.Rproj       
├── figures                          --  Contains outputs from simulations and performance comparison
│   ├── performance                  
│   └── simulations                  
├── renv                             --  [renv]
├── renv.lock                        --  [renv]
└── src
    └── estimate-degrees.cpp         --  Impute forward citations
```
