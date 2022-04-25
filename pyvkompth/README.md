# Instructions for Compilation of the Python Visual Interactive Wrapper for vKompth

## Requirements

0) You need a working `python3` installation, including `numpy` and `matplotlib`.

## Compilation

1) First compile the two model variants (bb, for blackbody and
  dk, for disk-blackbody) as python modules, editing `-lopenblas`
  modifier to link your desired library option, as follows:

```
python -m numpy.f2py -lopenblas -c pyvkompthbb.pyf pyvkompthbb.f90 \
       ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
       ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
       ../sco_model_LOGbb.f90 ../dependencies/*.f ../xsbbrd.f
```

```
python -m numpy.f2py -lopenblas -c pyvkompthdk.pyf pyvkompthdk.f90 \
       ../sco_arrays.f90 ../sco_global.f90 ../sco_band_integration.f90 \
       ../sco_mppinv.f90 ../sco_simpson.f90 ../sco_par.f90 \
       ../sco_model_LOG_dskb.f90 ../dependencies/*.f ../xsdskb.f
```

## Run the interactive wrapper

2) Before running the visually interactive python wrapper, please set your
`OPENBLAS_NUM_THREADS` variable to the amount of threads wanted.
Otherwise, it will use ALL the threads available.

3) Run the interactive wrapper as `python3 pyvkompth.py`
