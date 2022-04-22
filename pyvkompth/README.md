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
