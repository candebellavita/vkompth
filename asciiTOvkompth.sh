#!/bin/bash

if [ -z $HEADAS ]
then
	echo
	echo '  ERROR: Please, initiate HEADAS environment.'
	echo
	exit
fi

echo '   ======================================================='
echo '     This code translates ascii rms/lags into FITS files'
echo '     compatible with the variable-Comptonization model by'
echo '     Bellavita, Garcia, Mendez and Karpouzas (2022).'
echo '     Original model in: Karpouzas et al. (2020).'
echo '     Please cite these papers if you use this model.'
echo '     Feel free to contact us through email or GitHub.'
echo '   ======================================================='
echo '     Units:  rms in fractional units and lags in radians'
echo '     Format: E_lo E_hi [rms, lag] [rms_err,lag_err] (68%)'
echo '   ======================================================='
echo '     Usage: asciiTOvkompt file.ascii [rms, lag] QPOfreq.'
echo '   ======================================================='
echo

if [ $# -ne 3 ]
then
    echo
    echo '  ERROR: Please run with 3 arguments:'
    echo '              asciiTOvkompt file.ascii [rms, lag] QPOfreq.'
    echo
    exit
fi

ascii=$(basename $1 .ascii)

echo '    Converting to FITS...'
awk '{print $1,$2,$3*($2-$1),$4*($2-$1)}' ${ascii}.ascii > ${ascii}.$$.tmp
ftflx2xsp ${ascii}.$$.tmp ${ascii}.pha ${ascii}.rmf clobber=yes
rm ${ascii}.$$.tmp
echo '                         ... OK'

if [ $2 == 'rms' ]
then
	fparkey "mode: 1" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'

	echo ''
	echo '   Keywords Mode: 1 (rms) ; QPO: '$3' Hz'
        echo '              were succesfully added.'
        echo ''

elif [ $2 == 'lag' ]
then
	fparkey "mode: 2" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'

	echo ''
	echo '    Keywords Mode: 2 (lag) ; QPO: '$3' Hz'
        echo '              were succesfully added.'
        echo ''

elif [ $2 == 'sss' ]
then
	fparkey "mode: 3" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'

	echo ''
	echo '    Keywords Mode: 3 (sss) ; QPO: '$3' Hz'
        echo '              were succesfully added.'
        echo ''

elif [ $2 == 'real' ]
then
	fparkey "mode: 4" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'

	echo ''
	echo '    Keywords Mode: 4 (real) ; QPO: '$3' Hz'
        echo '              were succesfully added.'
        echo ''

elif [ $2 == 'imag' ]
then
	fparkey "mode: 5" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'

	echo ''
	echo '    Keywords Mode: 5 (imag) ; QPO: '$3' Hz'
        echo '              were succesfully added.'
        echo ''

else
	echo ''
	echo '  ERROR: Indicate spectrum type and QPO frequency.'
	echo ''
	exit
fi
