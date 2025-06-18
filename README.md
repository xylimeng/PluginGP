# PluginGP
Inference on function derivatives using Plug-in Gaussian processes

## Organization

### R/ (Source code)
- Contains source code for plug-in Gaussian process used in the paper.
  - `GP_Matern.R`: Plug-in GP with Matérn kernel.
  - `GP_SE.R`: Plug-in GP with squared-exponential kernel.
  - `GP_Sobolev.R`: Plug-in GP with Sobolev kernel.
  - `GP_Matern_Hetero.R`, `GP_SE_Hetero.R`, `GP_Sobolev_Hetero.R`: GP with heterogeneous error used for real data application Figure 2(b).
  - `GP_Bspline.R`: GP with B-spline prior in Yoo and Ghosal (2016) as a benchmark; the code was provided by the original author of the paper.

### simulation/ (Section 4 of the paper)
- #### code/
  - `Figure1.R`: Produces plots shown in Figure 1.
  - `Table1.R`: Generates results in Table 1.

- #### result/
  - Simulation outputs used to generate Table 1.

- #### figure/
  - Plots in Figure 1.

### real_data/ (Section 5 of the paper)
- #### code/
  - `Figure2a.R`: Code for generating Figure 2(a).
  - `Figure2b.R`: Code for generating Figure 2(b).

- #### data/
  - `CSIRO_Recons_gmsl_yr_2015.txt`: Global mean sea level data. The data used in our analysis are directly downloaded from \url{https://research.csiro.au/slrwavescoast/sea-level/measurements-and-data/sea-level-data}, under the section “Update of Reconstructed GMSL from 1880 to 2013” within “Sea level and ocean heat content.” The dataset is already cleaned and formatted for use; no additional preprocessing is required.


- #### figure/
  - `GMSL-1.pdf`: Corresponds to Figure 2(a).
  - `GMSL-2.pdf`: Corresponds to Figure 2(b).

## Dependencies

The following R packages are required to run the code:

- `splines`
- `MASS`
- `Rsolnp`
- `RandomFieldsUtils`
- `foreach`
- `doSNOW`
- `doParallel`

The code depends on the `RandomFieldsUtils` package for computing the Matérn kernel. Since this package has been removed from CRAN, it must be installed manually from GitHub:

```r
install.packages("remotes")  # if not already installed
remotes::install_github("cran/RandomFieldsUtils")
```