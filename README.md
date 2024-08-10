# vKompth

(It is "v-komp-t-h"; phonetically it is "viː-käm(p)-tiː-eɪʧ") 

This is the Fortran version of the [vKompth](https://github.com/candebellavita/vkompth) code from [Bellavita et al. 2022, MNRAS 515, 2099](https://ui.adsabs.harvard.edu/abs/2022MNRAS.515.2099B/) which was originally developed in [Karpouzas et al. 2020, MNRAS 492, 1399](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1399K/) and [García et al. 2021, MNRAS 501, 3173](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.3173G/).


## Requirements

In order to invert the matrix of the linear problem, `vKompth` requires the `DGTSV`, `ZGETRF` and `ZGETRS` routines from the [Lapack](https://www.netlib.org/lapack) or [openBLAS](https://www.openblas.net/) libraries, which have to be installed in your system. Before compiling the code check the `lopenblas` choice under the `LDLIBS` variable in the `Makefile` and change it to `lopenblas`, `lopenblasp` or `llapack`, depending the one you have installed in your system, and `LDFLAGS` to the directory where the library is installed.

Make sure you have [HEASOFT](https://heasarc.gsfc.nasa.gov/lheasoft/) installed from the source version, and set it up before compiling the *XSPEC* models. This model has to be compiled with the same compiler version used for *HEASOFT* to ensure compatibility. Read the *HEASOFT* pages about compilers.

In order to compile the *python* wrappers, you need to have a working `python3` environment set up, including `numpy` and `matplotlib` packages.


## Compile the main code, *XSPEC* models and *python* wrapper

Run `make` to compile the main program `vKompth` source code, together with the corresponding *XSPEC* models and *python* wrappers, including the different model variants: **bb=blackbody** seed-photon source; **dk=diskbb** seed-photon source; **dual=two coronas**.

To run `vKompth` in multithread mode set, for instance:
```
export OPENBLAS_NUM_THREADS=2  #(BASH version)
setenv OPENBLAS_NUM_THREADS 2  #(CSH version)
```


## Run the *python* wrapper

Go into `pyvkompth` subdirectory and run `python3 pyvkompth.py`. A GUI will load plotting both *rms* and *lags* for `vkompthbb` and `vkompthdk` model variants, allowing to modify parameters on-the-fly with interactive sliders.


## Load and run the *XSPEC* models

In an XSPEC session, for instance, `vkompthbb` can be then loaded using:
```
lmod vkompthbb /PATHTO/vkompthbb/
```
(and similar commands for the other wrappers: `vkompthdk`, `vkdualbb` and `vkdualdk`).

Alternatively, the four model variants can be loaded using the `XSPEC` script provided as `@load_vkompth.xcm`. This script also includes plotting commands like `plrl` which can be used to produce a fancy plot for both *rms* and *lags* in *XSPEC* if data are loaded as `data 1:1 rms.pha 2:2 lag.pha`.

Examples of those can be found under the `MAXI_J1348-630` subdirectory, which correspond to data published under [Bellavita et al. 2022, MNRAS 515, 2099](https://ui.adsabs.harvard.edu/abs/2022MNRAS.515.2099B/).

Under the directory `DILUTION` you can find an example of the simultaneous fit to the spectrum of the source and the rms and lag spectra of the QPO including the effect of dilution on the rms spectrum. See for instance [Ma et al. 2023, MNRAS, 525, 854](https://ui.adsabs.harvard.edu/abs/2023MNRAS.525..854M/), or [Zhang et al. 2023, MNRAS, 520, 5144](https://ui.adsabs.harvard.edu/abs/2023MNRAS.520.5144Z/) for details.


## Preparing your own data

We also provide a `bash` script named `asciiTOvkompth.sh`, which uses `FTOOLS` to convert ASCII files with energy-dependent *rms*, *lags* and *time-averaged* spectra, into *PHA* and *RMF* FITS files fully compatible with our `vKompth` model variants (by including both *QPO frequency* and *mode* as *XFLT* variables for *XSPEC*).


## Questions, comments, issues

For questions, comments and issues, please contact the authors of [Bellavita et al. 2022, MNRAS 515, 2099](https://ui.adsabs.harvard.edu/abs/2022MNRAS.515.2099B/) preferably through the [vKompth GitHub](https://github.com/candebellavita/vkompth).
