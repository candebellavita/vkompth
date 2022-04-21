# vKompth

This is the Fortran version of the [vKompth](https://github.com/candebellavita/vkompth) code from *Bellavita et al. 2022 (subm. to MNRAS)* which was originally developed in [Karpouzas et al. 2020, MNRAS 492, 1399](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1399K/abstract) and [García et al. 2021, MNRAS 501, 3173](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.3173G/abstract).


## Requirements

In order to invert the matrix of the linear problem, `vKompth` requires the `DGTSV`, `ZGETRF` and `ZGETRS` routines from the [Lapack](https://www.netlib.org/lapack) or [openBLAS](https://www.openblas.net/) libraries, which have to be installed in your system. Before compiling the code check the `lopenblas` choice under the `LDFLAGS` variable in the `Makefile` and change it to `lopenblas`, `lopenblasp` or `llapack`, depending the one you have installed in your system.

Make sure you have [HEASOFT](https://heasarc.gsfc.nasa.gov/lheasoft/) installed from the source version, and set it up before compiling the XSPEC models. This model has to be compiled with the same compiler version used for HEASOFT to ensure compatibility. Read the HEASOFT pages about compilers.


## Compile the code

Run `make` to compile the main program `vKompth` source code.

Alternatively, to manually compile the code, use, for instance:
```
gfortran -O5 -Wall dependencies/*f sco_simpson.f90 sco_mppinv.f90 sco_model.f90 sco_band_integration.f90 sco_par.f90 sco_arrays.f90 sco_global.f90 sco_programDSKB.f90 -lopenblas -o vkompth_dk
```
replacing `lopenblas` by the appropriate library for your system.

To run `vKompth` in multithread mode use, for instance:
```
export OPENBLAS_NUM_THREADS=2  #(BASH version)
setenv OPENBLAS_NUM_THREADS 2  #(CSH version)
```


## XSPEC models or wrappers

You will now need to compile the different versions of the program (the so-called wrappers); the names of the subdirectories indicate the type of wrapper model: **bb=blackbody** seed-photon source; **dk=diskbb** seed-photon source; **dual=two coronas**.

To compile the XSPEC wrappers, go to each wrapper subdirectory and, after loading HEASOFT, run `initpackage` and `hmake`. For instance:
```
cd vkompthbb
```
Please, check and manually edit `Makefile_libs` to match the corresponding `lopenblas` value for your system if necessary. Then run:

```
initpackage vkompthbb lmod_vkompthbb.dat  . ; cp Makefile_libs Makefile; hmake
```

Follow similar steps for the other wrappers.

In an XSPEC session, `vkompthbb` can be then loaded using:
```
lmod vkompthbb /PATHTO/vkompthbb/
```
(and similar commands for the other wrappers).


## Questions, comments, issues

For questions, comments and issues, please contact the authors of *Bellavita et al. 2022 (subm. to MNRAS)* preferably through the [vKompth GitHub](https://github.com/candebellavita/vkompth).
