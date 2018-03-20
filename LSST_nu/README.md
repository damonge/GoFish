## LSST/CMB-S4 complementarity for neutrino masses and dark energy

This folder contains the files necessary to reproduce the results of 1803.xxxxx studying the complementarity between CMB-S4 and LSST to distinguish between neutrino masses and non-minimal cosmological scenarios. The following files are provided:
- `param_LSST_DESI_Cls.ini` and `PlotCls.ipynb`: GoFish configuration file and Jupyter notebook to reproduce the plots showing binned power spectra (Fig. 2).
- `param_LSST_DESI_Nu.ini` and `PlotTrianglesTables.ipynb`: GoFish configuration file and Jupyter notebook to reproduce all other plots in the paper.

The `plots` folder contains all output figures, and `run_class.sh` runs the spectra in parallel as SLURM jobs on the Princeton clusters if specified in `exec_path` in the config files.  