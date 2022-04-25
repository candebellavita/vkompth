# vKompth

This is the Fortran version of the [vKompth](https://github.com/candebellavita/vkompth) code from *Bellavita et al. 2022 (subm. to MNRAS)* which was originally developed in [Karpouzas et al. 2020, MNRAS 492, 1399](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1399K/abstract) and [García et al. 2021, MNRAS 501, 3173](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.3173G/abstract).


## Requirements

In order to invert the matrix of the linear problem, `vKompth` requires the `DGTSV`, `ZGETRF` and `ZGETRS` routines from the [Lapack](https://www.netlib.org/lapack) or [openBLAS](https://www.openblas.net/) libraries, which have to be installed in your system. Before compiling the code check the `lopenblas` choice under the `LDFLAGS` variable in the `Makefile` and change it to `lopenblas`, `lopenblasp` or `llapack`, depending the one you have installed in your system.

Make sure you have [HEASOFT](https://heasarc.gsfc.nasa.gov/lheasoft/) installed from the source version, and set it up before compiling the *XSPEC* models. This model has to be compiled with the same compiler version used for *HEASOFT* to ensure compatibility. Read the *HEASOFT* pages about compilers.

In order to compile the *python* wrappers, you need to have a working `python3` environment set up, including `numpy` and `matplotlib` libraries.


## Compile the main code, *XSPEC* models and *python* wrapper

Run `make` to compile the main program `vKompth` source code, together with the corresponding *XSPEC* models and *python* wrappers, including the different model variants: **bb=blackbody** seed-photon source; **dk=diskbb** seed-photon source; **dual=two coronas**.

To run `vKompth` in multithread mode set, for instance:
```
export OPENBLAS_NUM_THREADS=2  #(BASH version)
setenv OPENBLAS_NUM_THREADS 2  #(CSH version)
```


## Load and run the *XSPEC* models

In an XSPEC session, for instance, `vkompthbb` can be then loaded using:
```
lmod vkompthbb /PATHTO/vkompthbb/
```
(and similar commands for the other wrappers: `vkompthdk`, `vkdualbb` and `vkdualdk`).


## Run the *python* wrapper

Go into `pyvkompth` subdirectory and run `python3 pyvkompth.py`. A GUI will load plotting both *rms* and *lags* for `vkompthbb` and `vkompthdk` model variants, allowing to modify parameters on-the-fly with intereactive sliders.


## Questions, comments, issues

For questions, comments and issues, please contact the authors of *Bellavita et al. 2022 (subm. to MNRAS)* preferably through the [vKompth GitHub](https://github.com/candebellavita/vkompth).
