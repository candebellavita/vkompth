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
