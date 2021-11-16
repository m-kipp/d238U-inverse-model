# d238U inverse model

Version 1.0

This `R` code executes the inverse analysis of d238U paleo-redox datasets as described in Kipp & Tissot (2022) EPSL. 

Link to paper: https://doi.org/10.1016/j.epsl.2021.117240 

As a matter of courtesy, we request that those using this code please cite Kipp & Tissot (2022). In the interest of an "open source" approach, we also request that authors who use and modify the code please send a copy of papers and modified code to the lead author (mkipp@caltech.edu).

REQUIREMENTS: `R` (this code was written using v4.0.0), including `msir`, `FME`, `doParallel` and `LaplacesDemon` packages (and their dependencies).

HOW TO RUN CODE:
1) Open the `R` script (both the `.R` and `.Rmd` files work, and both can be opened in `R` or `R Studio`).
2) Set working directory, and read a d238U dataset from that working directory. *Nota bene*: ensure that the column names are 'time', 'd238U' and 'err' ('err' should be the 1SD analytical uncertainty).
3) Select model parameters (`time_step`, `prop_uncert`, `m`, `niterMCMC`, `updatecovMCMC`, `n_walkers`).
4) Execute the entire code. One can do this piece-by-piece and check for errors, or can run the full analysis and assess the outcome afterward. The final retrieved trends will be plotted at the end, along with relevant statistics for assessing convergence of the MCMC routine. If the model has not converged, the model parameters can be adjusted and the code executed again, iterating until it converges upon the best-fit solution.

An annotated walk-through of this workflow is provided in the attached R Markdown document (`d238U-inverse-model.Rmd`) and viewable at https://m-kipp.github.io/d238U-inverse-model/R-Markdown. The dataset used in the walk-through (`jost_2017_Triassic_data.csv`) is also available for download, to allow a test run before trying the model on a new dataset. 
