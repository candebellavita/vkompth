# Instructions for Compilation under XSPEC

1) First load HEASOFT (heainit, export, etc...).
2) Edit the Makefile_libs LDFLAGS and LDLIBS variables
   to your preferred values (lopenblas, lopenblasp, llapack, etc).
3) Then, run:

```
initpackage vkdualdk lmod_vkdualdk.dat .
cp Makefile_libs Makefile
hmake
```

We use Makefile_libs to inject all the dependencies of our
XSPEC model with the files on the top dir, and openblas
library for multithreading.

3) Before running XSPEC, please set your `OPENBLAS_NUM_THREADS`
variable to the amount of threads wanted. Otherwise, it will
use ALL the threads available, which may interfere with your
parallel sets under XSPEC.

4) Please edit and use `load_vkompth.xcm` to load the model into XSPEC together with some hacks to plot fractional rms, and phase lags.

5) The model can also be compiled as a python module using .pyf.

```
python -m numpy.f2py -lopenblas -c vkdualdk.pyf vkdualdk.f90 \
       ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
       ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
       ../sco_model_LOG_dskb.f90 ../dependencies/*.f ../xsdskb.f \
       -L$HEADAS/../Xspec/BLD/x86_64-pc-linux-gnu-libc2.31/lib/ -lXSUtil -lXSFunctions
```
