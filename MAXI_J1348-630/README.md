# Example of how to use the model

We assume that you followed the instructions and you were able to compile and load the model.

We provide a number of files under the subdirectory `MAXI_J1348-630/` that you can use to experiment with the model. We first show below how to load the data, define and fit the model, and plot the results. We first fit a 1-corona and afterwards a 2-corona model to the data. At the end we explain how to produce rms and lag spectra that are suitable to fit with the model.

## 1. Example of fitting with 1 corona:

```
XSPEC12>lmod vkompthdk /full_path_to_the_directory_with_the_model/
```

Use the .xcm file `@vkompthdk_init.xcm` to load the data (data 1:1 is the rms and data 2:2 is the lags), define the model and give some initial values to the parameters:

```
XSPEC12>@vkompthdk_init.xcm
…
XSPEC12>show fit

Fit statistic  : Chi-Squared                43051.45     using 16 bins.
                 Chi-Squared                 1786.84     using 16 bins.
Total fit statistic                         44838.29     with 27 d.o.f.

Test statistic : Chi-Squared                44838.29     using 32 bins.
 Null hypothesis probability of 0.00e+00 with 27 degrees of freedom
 Current data and model not fit yet.
```

Notice that the parameters for the second dataset (the lags in this case) are the same as (linked to) the ones of the first dataset (the rms in this case); the data files will contain a keyword (see below) that tells the model whether it has to output the rms or lags vs. energy.

Because the rms and lag are relative quantities, the model cannot be renormalised. It is therefore a good idea to plot the data and the fit and make some small adjustment to 2 of the parameters to make sure that the model on average matches the data. We make use of the new functionality of XSPEC in HEASOFT 6.30 to plot multiple datasets separately

```
XSPEC12>set plot ene; cpd /xs
XSPEC12>plot 1 lda de 2 da de

```

We see that there is and offset between the model both for the rms and the lags. Show the free parameters to see which ones to change:

```
XSPEC12>show fre

Free parameters defined:
========================================================================
Model vkompthdk<1> Source No.: 1   Active/On
Model Model Component  Parameter  Unit     Value
 par  comp
                           Data group: 1
   1    1   vkompthdk  kTs        keV      0.500000     +/-  0.0          
   4    1   vkompthdk  size       km       1000.00      +/-  0.0          
   5    1   vkompthdk  eta                 0.500000     +/-  0.0          
   8    1   vkompthdk  DHext               5.00000E-02  +/-  0.0          
  10    1   vkompthdk  reflag              0.0          +/-  0.0          
                           Data group: 2
________________________________________________________________________

```

To make the rms of the model roughly match that of the data reduce the parameter `dHext` by a factor of ~10; this will reduce the average rms in the model by the same factor:

```
XSPEC12>new 8 .05
```

You can replot to confirm.

Next adjust reflag such that the average lag of the model is more or less the same as that of the data. In this case the lags of the model are too low compared to the data, so increase `reflag` by ~0.2

```
XSPEC12>new 10 .2
```

Again, replete to confirm.

You can now fit. You should get:

```
XSPEC12>show fre

Free parameters defined:
========================================================================
Model vkompthdk<1> Source No.: 1   Active/On
Model Model Component  Parameter  Unit     Value
 par  comp
                           Data group: 1
   1    1   vkompthdk  kTs        keV      0.444977     +/-  3.35827E-02  
   4    1   vkompthdk  size       km       1.12952E+04  +/-  3279.69      
   5    1   vkompthdk  eta                 0.361139     +/-  6.45261E-02  
   8    1   vkompthdk  DHext               0.130668     +/-  3.58378E-02  
  10    1   vkompthdk  reflag              0.122488     +/-  9.19909E-03  
                           Data group: 2
________________________________________________________________________
```

and 

```
XSPEC12>show fit

Fit statistic  : Chi-Squared                   12.36     using 16 bins.
                 Chi-Squared                   80.80     using 16 bins.
Total fit statistic                            93.16     with 27 d.o.f.

Test statistic : Chi-Squared                   93.16     using 32 bins.
 Null hypothesis probability of 3.32e-09 with 27 degrees of freedom
```

Plot and see the result.

We also provide an `.xcm` file with our best-fitting model. You can load it and compare with your result:

```
XSPEC12>@vkompthdk.xcm 
…
XSPEC12>show free
…
XSPEC12>show fit
…
XSPEC12>plot
… etc.
```

While the fit looks relatively okay, as in Garcia et al. (2021), the errors of the rms have been multiplied by a factor so that the rms does not (totally) dominate the fit. Besides, there are some systematic negative residuals both in the rms and the lags at around 2-3 keV, and positive residuals at low and high energies.

Let’s therefore try a 2-corona model.

## 2. Example of fitting with 2 coronas: 

Given that the models with 2 coronas (dual models) have more free parameters than those of a single corona, you need to have enough data points (rms plus lags) compared to the number of parameters.

```
XSPEC12>lmod vkdualdk /full_path_to_the_directory_with_the_model/
```

Use the .xcm file `@vkdualdk_init.xcm` to load the data (data 1:1 is the rms and data 2:2 is the lags), define the model and give some initial values to the parameters:

```
XSPEC12>@vkdualdk_initial.xcm 
```

To start the fit of the 2 coronas at the point where we ended with a 1-corona model, in this `.xcm` file we set all the parameters of the second corona equal to those of the first corona and the phase of the oscillations of the second corona with respect to the first one equal to 0. By doing this, the dual-corona model reproduces the case of a single corona. (Notice that in this case the rms values are not multiplied by the factor but have the true errors, and hence the chi^2 is bigger than in the previous example.)

We know untie the parameters of the second corona (you can check the number of the parameters with `show all` or `show tied`) and free the phase between the 2 coronas:

```
XSPEC12>untie 2 4 6 8 10 14 ; thaw 15
```

It is now wise to give a small “kick” to the parameters of the coronas to get the model out of the current minimum of the chi^2 (if you do not do that, it may stay there and miss a deeper minimum). E.g.:


```
XSPEC12>new 1 .6; new 7 1000; new 9 .8; new 9 .1; new 15 1
```

Plot the data and model and change `dHext1`, `dHext2` and `reflag` to match the model to the data on average. In this case it is a bit less intuitive how to do that since the effect of the two `dHext` parameters depends also on the value of `phi`, the relative phase of the 2 coronas. For instance, in this case I get a reasonable match if I set `dHext1` to 0.25 and do not change `reflag` (you will have to experiment in other cases):

```
XSPEC12>new 13 .25
```

Then fit and after a while you should get:

```

XSPEC12>show fre

Free parameters defined:
========================================================================
Model vkdualdk<1> Source No.: 1   Active/On
Model Model Component  Parameter  Unit     Value
 par  comp
                           Data group: 1
   1    1   vkdualdk   kTs1       keV      0.958184     +/-  0.501929     
   2    1   vkdualdk   kTs2       keV      0.423242     +/-  8.94960E-02  
   7    1   vkdualdk   size1      km       163.448      +/-  430.400      
   8    1   vkdualdk   size2      km       1.18282E+04  +/-  9200.07      
   9    1   vkdualdk   eta1                0.965229     +/-  0.180595     
  10    1   vkdualdk   eta2                0.214920     +/-  0.650714     
  13    1   vkdualdk   DHext1              0.204780     +/-  0.237245     
  14    1   vkdualdk   DHext2              0.211659     +/-  0.146753     
  15    1   vkdualdk   phi                 -3.26031     +/-  0.384806     
  17    1   vkdualdk   reflag              0.145333     +/-  1.31813E-02  
                           Data group: 2
________________________________________________________________________

XSPEC12>show fit

Fit statistic  : Chi-Squared                   13.40     using 16 bins.
                 Chi-Squared                    7.68     using 16 bins.
Total fit statistic                            21.08     with 22 d.o.f.

Test statistic : Chi-Squared                   21.08     using 32 bins.
 Null hypothesis probability of 5.16e-01 with 22 degrees of freedom

```

Plot to see the fit.

You can compare your fit with the one we got (use file `@vkdualdk.xcm`).

As usual, in other cases you may have to run `error` and `steppar` a few times to see if you can find a deeper minimum. It is also a good idea to run a long MCMC to see if the fit converges to a different solution, or whether the posterior of the parameters is multimodal.

## 3. Making the rms and lag spectra necessary for the fits:

You will need .pha/.rmf pairs for the rms and the lag spectra. You can make them with your own tools, or use the tool provided here. If you use your own tool, you will need to add two keywords to the rms and lag spectra before you fit. We explain this at the end. If you use the tool provided here to make the files the keywords will be added automatically.

The rms spectrum must be in fractional units, with values between 0 and 1.

The lag spectrum must be phase lags in radians, with values from -pi/2 to pi/2.

Since the lags are a relative measurement (lags of a subject band with respect to some reference band), the reference band for the lags is not important because this only adds a constant to the lag spectrum; the model accounts for this (this is equivalent to the normalisation in additive models). You can therefore choose any band as the reference band. Notice that if you use the full band as reference, and this one encompasses some of the subject bands, you will need to correct the lags for the correlation introduced by the photons that are the same in both the subject and reference band. See Ingram 2019 for that.

If you use any other band as reference, for instance, the lowest energy band, and you want to include that band (with 0 lags) in the fits, you need to assign an error to that lag. For instance, you can use the average error of all the other bands. 

To create the rms and lag spectra using the tool provided, you will need ASCII files with the data (units explained above) that you will use to make the pha+rmf files (see below). These ASCII files must have an extension ‘.ascii’ in the name; they could be, for instance, ‘rms_data.ascii’ and ‘lag_data.ascii’.

The ASCII file for the rms spectrum of the QPO should have 4 columns:

Emin Emax fractional_rms 1-sigma_error
…

where Emin and Emax are the minimum and maximum energy of the band for which you have measured the fractional rms amplitudes (with errors), and each row is a new measurement.

The ASCII file for the phase-lag spectrum of the QPO should have 4 columns:

Emin Emax phase-lag_in_rad 1-sigma_error
…

where Emin and Emax are the minimum and maximum energy of the band for which you have measured the lags (with errors), and each row is a new measurement.

The bands for the rms and the lags do not have to be the same.

To make the .pha and .rmf files use the command:

```
asciiTOvkompth rms_data rms QPO_frequency_in_Hz
asciiTOvkompth lag_data lag QPO_frequency_in_Hz
```

Notice that you do not have give the extension of the ASCII files, which is assumed to be ‘.ascii’.

In both cases we give the frequency of the QPO in Hz. For example, assuming that the QPO frequency is 4.5 Hz, the commands:

```
asciiTOvkompth rms_data rms 4.5
asciiTOvkompth lag_data lag 4.5
```

take the files `rms_data.ascii` and `lag_data.ascii` and create the files `rms_data.pha / rms_data.rmf’ and `lag_data.pha / lag_data.rmf`, respectively.

The .pha files will have a keyword in the header indicating whether the data are rms amplitudes or lags, and another keyword with the frequency of the QPO given in the command line. These keywords will be used by the model to compute the model of either the rms or the lag for the given frequency.

## 4. Notes:

Do not combine 2 or more instances of this model in Xspec. The correct way to do this is to add components in the Fourier, Real and Imaginary, space. Adding components may improve the fits, but will have no physical meaning.

Do not add other additive components (body, powerlaw, etc .) to the model. The same considerations mentioned in the previous point apply here.

Do not add multiplicative components (phabs, gabs, etc.) to the model. The same considerations mentioned in the previous point apply here.

Notice, by the way, that absorption does not affect the rms or lags.

Although the rms amplitude is independent of the effective area of the instrument that is used to measure them, the background does affect the rms. This is so because the rms is defined in terms of the power, P, the total observed count rate, C, and the background count rate, B, as:

`rms =\sqrt{P/C} * (C)/(C - B)`.

If the total and background count rates used in the calculation are wrong (e.g., the background comes from a model but is not properly calibrated, the total count rate includes contamination from other sources that is not taken into account in the background count rate, or there are particle flares during the observation that are not accounted for), the rms may be biased one way or another. This may be noticeable if you are fitting simultaneous observations with different instruments. (Notice that the same may happen if you measure the same source with the same instruments at different times. The rms amplitude of the QPO of the source may not change, but a background flare, not properly accounted for, may affect the observed rms amplitude.)


## 5. More advanced usage of the model:

### 5.1 Fitting rms and lags of the QPO and source spectrum simultaneously.

You can use these models to fit also the rms and lags of the QPO and the total energy spectrum simultaneously with common parameters across model components linked. In this case you read 3 datasets, e.g.:

```
data 1:1 rms.pha 2:2 lag.pha 3:3 source.pha
```

You can then define a model like:

```
model vkompthdk + phabs*(diskbb+vkompthdk)
```

As before, the rms and lags will have a keyword that is passed to the model; since the source spectrum does not have such keyword the model `vkompthdk` will output the time-averaged total spectrum (simular to `nthcomp`).

One has to be careful when defining the parameters such that the first `vkompthdk` component fits the rms and lags of the QPO and the rest of the model fits the source spectrum. You can then link the common parameters between the first and second `vkompthdk`and the temperatures of the seed photon source in `diskbb` and the two `vkompthdk`. Since the time-averaged version of `vkompthdk` is equivalent to `nthcomp`, you can use the latter if you prefer instead of the second `vkompthdk`.

One should then fix `NH` of `phabs` and the norms of `diskbb` and the second `vkompthdk` (or `nthcomp`) to 0 for the data sets 1 and 2  (rms and lag spectra) and the norm of the first `vkompthdk` to 0 for the dataset 3 (time-average source spectrum). 

It is possible to simplify the model further by using the second `vkompthdk` component both for the rms and lags on one hand and the total spectrum on the other. In that case you will need to be careful on how you define the parameters to make sure that you fit the right model to the right data.

We leave this explanation here. Please contact us if you have trouble implementing this case. 

### 5.2 Dilution.

If, as this model assumes, the variability in the QPO data comes only from the corona while all other additive components in the spectrum are not variable, those other components will dilute the variability of the corona in the observed rms spectrum as a function of energy. For instance, if there is a non-variable disc in the spectrum that dominates the emission below, say, 2 keV, the observed rms amplitude below 2 keV will be lower than the intrinsic rms amplitude of the corona (which is what the model computes) at those energies. Because of this, one needs to correct the rms amplitude given by the model for this dilution effect when one compares the data with the model. (N.B.: dilution does not affect the lags.)

If you obtained a fit of the total spectrum separately, you can apply a dilution factor during the fits, fixing the parameters of the components in the energy spectrum to your best-fitting values to calculate the dilution.  Alternatively (and possibly more correctly), you can fit the rms, lags and total spectrum simultaneously and calculate the dilution on the fly during the fits. 

A quick (and not complete) explanation of how to do this is:

a. Read your data and define your model as in 5.1.
b. Use the command `mdefine` in Xspec to define the dilution. For example:

`XSPEC12>mdefine dilution nthcomp(Tin,Gam,kTe,1,0)*n_nth / (diksbb(Tin)*n_dbb + nthcomp(Tin,Gam,kTe,1,0)*n_nth) : mul`

where we assumed the disc version of `nthcomp` with 0 redshift; notice also that you cannot use `Gamma` for the name of a parameter because Xspec has the `gamma` function and this would confuse the program.

(We did not include `phabs` in the `dilution` because multiplicative components cancel out in the ratio and do not contribute to the dilution.)

You can now add the component `dilution` to the model:

`XSPEC12> editmod vkompthdk*dilution + phabs*(diskbb+vkompthdk)`

and link the parameters as necessary.

We do not explain this in detail here. If you need help, please contact us.
