Version 1.2.1
=============

Added `PDBVAL` variables, which allow to interact through
`tcloutr modkeyval` to obtain `eta_int` and their uncertainties
using `simpars`.

Added `etaint.tcl` and `etaint12.tcl` TCL XSPEC scripts to
calculate `eta_int` and their uncertainties using `simpars`.


Version 1.2.0
=============

Added `dL` variant including perturbations to the size of the corona,
electron density and optical depth, which replaces the variation in
external-heating rate.

Fixed bug in the calculation of `dkTs` and `dkTe` internal variables.
Values in previous versions are unreliable. 

Fixed insignificant bug in the application of the Klein Nishina correction.

`phi` is now defined as periodic in `(-pi, pi)` in the `lmod_.dat` of dual variants.


Version 1.1.2
=============

Added `eta_int` and `outflux` variables to `xset`
within the XSPEC wrappers.

Changed `Emin_adim` to 1e-6 and meshlog to 494 to
keep numerical resolution.


Version 1.1.1
=============

Added `xset` *eta* variables to be able to recover the
internal value of `\tilde\eta` for each model version.

Added `power` mode, to match Lorentzian normalizations
fitted with XSPEC to the QPOs in the PDS.


Version 1.1.0
=============

Updated vKompth XSPEC wrappers to include QPO frequency
and Mode (rms, lag, sss, etc) as XFLT internal variables.

For this purpose, `asciiTOvkompth.sh` was added to help
creatings FITS PHA and RMF files from ASCII files containing
rms and lag data.

Note that each version of the model now has 2 parameters less,
thus `xcm` model files are not backwards/forward compatible.

vKompth now prints a message indicating the model version
when first loaded in an XSPEC sesion.

Added a Python Interactive-Plot wrapper: `pyvkompth`

Makefile is now able to fully compile both VKOMPTH model,
python interactive-plot wrapper, and XSPEC wrappers when
HEASOFT is loaded.


Version 1.0.1
=============

Updated guidelines in `README.md` files for installation.

Updated Makefile code names.

Added `load_vkompth.xcm` XSPEC script to ease loading the
models and making quick pre-formatted plots.


Version 1.0.0
=============

Working versions of `vkompth[bb,dk]` and `vkdual[bb,dk]`.
- Use MODE parameter to select between rms, lag or sss spectrum.
- Use QPO parameter to input the QPO frequency in Hz.

MAXI J1348-630 includes pha/rmf data, and best-fitting xcm files.

Follow guidelines in `README.md` file for installation.
