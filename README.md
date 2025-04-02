# PluginGP
Inference on function derivatives using Plug-in Gaussian processes

## Organization

### R/ (Source code)
- Contains source code for plug-in Gaussian process used in the paper.
  - `GP_Matern.R`: Plug-in GP with Matérn kernel.
  - `GP_SE.R`: Plug-in GP with squared-exponential kernel.
  - `GP_Sobolev.R`: Plug-in GP with Sobolev kernel.
  - `GP_Bspline.R`: GP with B-spline prior in Yoo and Ghosal (2016).
  - `GP_Matern_Hetero.R`, `GP_SE_Hetero.R`, `GP_Sobolev_Hetero.R`: GP with heterogeneous error used for real data application Figure 2(b).

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
  - `CSIRO_Recons_gmsl_yr_2015.txt`: Global mean sea level data.

- #### figure/
  - `GMSL-1.pdf`: Corresponds to Figure 2(a).
  - `GMSL-2.pdf`: Corresponds to Figure 2(b).
