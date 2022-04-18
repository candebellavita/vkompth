# vKompth

This is the Fortran version of the [vKompth](https://github.com/candebellavita/vkompth) code from
*Bellavita et al. 2022 (subm. to MNRAS)* which was originally
developed in [Karpouzas et al. 2020, MNRAS 492, 1399](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1399K/abstract) and [García et al. 2021, MNRAS 501, 3173](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.3173G/abstract).

In order to invert the matrix of the linear problem, `vKompth` requires `dgtsv` routine from the [Lapack](https://www.netlib.org/lapack) or [openBLAS]([https://www.openblas.net/) libraries which have to be installed in your system. Before compiling the code, first check `lopenblas` choice under `LDFLAGS` variable in the `Makefile` (and change it to `lopenblas`, `lopenblasp`, `llapack`, according to your system).

Then, run `make` to compile the `vKompth` source code.

Alternative, to manually compile the code, use:
```
gfortran -O5 -Wall dependencies/*f sco_simpson.f90 sco_mppinv.f90 sco_model.f90 sco_band_integration.f90 sco_par.f90 sco_arrays.f90 sco_global.f90 sco_program.f90 -lopenblas -o scorpio_fortran
```
replacing `lopenblas` by the appropriate library for your system.

To run `vKompth` in multithread mode, use, for instance:
```
export OPENBLAS_NUM_THREADS=2
```

To compile the XSPEC wrappers, enter to each wrapper subdirectory and, after loading HEASOFT, run `initpackage` and `hmake` as follows:
```
cd vkompthbb
initpackage vkompthbb lmod_vkompthbb.dat  . ; cp Makefile_libs Makefile; hmake
```
Please, first check and manually edit `Makefile_libs` to match the corresponding `lopenblas` value for your system.

In an XSPEC session, `vkompthbb` can be then loaded using:
```
lmod vkompthbb /PATHTO/vkompthbb/
```

For questions, comments and issues, please contact the
authors of *Bellavita et al. 2022 (subm. to MNRAS)* preferably through the [vKompth GitHub](https://github.com/candebellavita/vkompth).
