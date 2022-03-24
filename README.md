# vKompth

This is the Fortran version of the vKompth code from
Bellavita et al. 2022 (submitted) which was originally
developed in Karpouzas et al. 2020.

To compile, first check the lopenblas FLAG in the Makefile
(and change it to lopenblas, lopenblasp, etc, as in your system).

Then, run:

```
make clean
make
```

Alternative, to compile, use (after checking for lopenblas, lopenblasp, etc):
```
gfortran -O5 -Wall dependencies/*f sco_simpson.f90 sco_mppinv.f90 sco_model.f90 sco_band_integration.f90 sco_par.f90 sco_arrays.f90 sco_global.f90 sco_program.f90 -lopenblas -o scorpio_fortran
```

To run in multithread mode, use, for instance:
```
export OPENBLAS_NUM_THREADS=2
```

To compile the XSPEC wrappers, use, for instance:
```
cd vkompthbb
hmake clean; initpackage vkompthbb lmod_vkompthbb.dat  . ; cp Makefile_libs Makefile; hmake
```
(please, first check and edit Makefile_libs with the corresponding lopenblas value)

Then, vkompthbb can be loaded into XSPEC using:
```
lmod vkompthbb /PATHTO/vkompthbb/
```

For questions, comments and issues, please contact the
authors of Bellavita et al. 2022 (subm. to MNRAS).
