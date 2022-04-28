#!/bin/bash

if [ -z $HEADAS ]
then
	echo
	echo '  ERROR: Please, initiate HEADAS environment.'
	echo
	exit
fi

echo '   ========================================================='
echo '     This code converts ascii rms/lags into FITS files'
echo '     compatible with the variable-Comptonization model by'
echo '     Bellavita, Garcia, Mendez and Karpouzas (2022).'
echo '     Original models: Karpouzas+2020, Garcia+2021.'
echo '     Please cite these papers if you use this model.'
echo '     Feel free to contact us through email or GitHub.'
echo '   ========================================================='
echo '     Units:  rms in fractional units [0,1]'
echo '             lags in radians [-pi/2,pi/2]'
echo '             QPO freq in Hz'
echo '     Format: E_lo E_hi [rms, lag] [rms_err,lag_err] (68%)'
echo '   ========================================================='
echo '     Usage: '
echo '       asciiTOvkompth filename_ascii [rms, lag] QPOfreq'
echo '     Examples: '
echo '       asciiTOvkompth some_file.dat    rms 4.5'
echo '       asciiTOvkompth other_file.ascii lag 3.1'
echo '   ========================================================='
echo

if [ $# -ne 3 ]
then
    echo
    echo '  ERROR: Please run with 3 arguments:'
    echo '    asciiTOvkompth filename_ascii [rms, lag] QPOfreq'
    echo
    exit
fi

file=$1

if [ ! -f ${file} ]; then
	  echo
    echo '  ERROR: File "'${file}'" not found!'
		echo
    exit
fi


ext=${file##*.}                 # file extension
ascii=$(basename $file .$ext)   # file basename

echo '    Converting to FITS...'
awk '{print $1,$2,$3*($2-$1),$4*($2-$1)}' ${file} > ${ascii}.$$.tmp
ftflx2xsp ${ascii}.$$.tmp ${ascii}.pha ${ascii}.rmf clobber=yes
rm ${ascii}.$$.tmp
echo '                         ... OK'

if [ $2 == 'rms' ]
then
	fparkey "mode: 1" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
        if [ $? -ne 0 ]; then
					 echo
           echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi
	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'
        if [ $? -ne 0 ]; then
					 echo
           echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi

	echo ''
	echo '   Keywords Mode: 1 (rms) ; QPO: '$3' Hz'
  echo '              were succesfully added.'
  echo ''

elif [ $2 == 'lag' ]
then
	fparkey "mode: 2" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
        if [ $? -ne 0 ]; then
					 echo
           echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi

	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'
        if [ $? -ne 0 ]; then
					 echo
           echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi


	echo ''
	echo '    Keywords Mode: 2 (lag) ; QPO: '$3' Hz'
  echo '              were succesfully added.'
  echo ''

elif [ $2 == 'sss' ]
then
	fparkey "mode: 3" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
        if [ $? -ne 0 ]; then
           echo
					 echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi

	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'
        if [ $? -ne 0 ]; then
           echo
					 echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi


	echo ''
	echo '    Keywords Mode: 3 (sss) ; QPO: '$3' Hz'
  echo '              were succesfully added.'
  echo ''

elif [ $2 == 'real' ]
then
	fparkey "mode: 4" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
        if [ $? -ne 0 ]; then
					 echo
           echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi

	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'
        if [ $? -ne 0 ]; then
					 echo
           echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi


	echo ''
	echo '    Keywords Mode: 4 (real) ; QPO: '$3' Hz'
  echo '              were succesfully added.'
  echo ''

elif [ $2 == 'imag' ]
then
	fparkey "mode: 5" ${ascii}.pha XFLT0001 add=yes comm='Spectrum type'
        if [ $? -ne 0 ]; then
           echo
					 echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi

	fparkey "QPO: "$3"" ${ascii}.pha XFLT0002 add=yes comm='QPO frequency [Hz]'
        if [ $? -ne 0 ]; then
					 echo
           echo '  ERROR: fparkey command failed on file' ${ascii}.pha
					 echo
           exit
        fi


	echo ''
	echo '    Keywords Mode: 5 (imag) ; QPO: '$3' Hz'
  echo '              were succesfully added.'
  echo ''

else
	echo ''
	echo '  ASCII file was converted to pha/rmf but KEYWORDS were not added to the spectrum.'
  echo '  The vKompth XSPEC model will not be able to identify the type of data.'
	echo '  ERROR: Indicate spectrum type and QPO frequency.'
	echo ''
	exit
fi
