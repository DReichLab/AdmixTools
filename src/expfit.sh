#!/bin/sh

parfile=$1
BIN=/home/np29/biology/Admix/src/../bin

if [ $# -lt 1 ]; then
  echo -e "USAGE: expfit.sh parfile"
  exit 1
fi
BOOL_T=1
BOOL_F=0

function verify_file
{
  if [ ! -s $1 ]; then
    echo "ERROR: $1 does not exist or is empty"
    return $BOOL_F
  fi
  return $BOOL_T
}

verify_file $parfile
if [ $? == $BOOL_F ]; then exit 1; fi


echo -e "Running rexpfit.r "
input=`cat ${parfile} | grep output: | awk '{print $2}'`
echo -e "Input file: " ${input}

output="expfit_${input}"
echo -e "Output file: " ${output}

# Set defaults
output_col=4   		# data column containing weighted correlation
lval=0.5		# start fit at lval genetic distance. Distances <lval ignored
hval=100		# end fit at hval genetic distance. Distances > hval ignored                	
affine=TRUE	 	# affine = TRUE/FALSE
plot=FALSE    	        # plot the output = TRUE/FALSE
jackknife=FALSE 	#jackknife = TRUE/FALSE

# specified in parfile
jackknife=`cat ${parfile} | grep jackknife: | awk '{print $2}'`

# Fit exponential using R 
echo -e "Expfit logfile:" ${output}".flog"
plot="TRUE"
echo -e "Output plot:" ${output}".pdf"
Rscript $BIN/rexpfit.r ${input} ${output} ${output_col} ${lval} ${hval} ${affine} ${plot} "FALSE"

if [ "${jackknife}" != 'NO' ]; then
	# one chromosome removed
	jOutput="expfit_${input}"
	echo -e "Jackknife logfile:" ${jOutput}".log"
	echo -e "Jackknife Results: \n" > ${jOutput}.log
	plot="FALSE"

	for chr in {1..22}
	do
		inputChr="${input}:${chr}"
		Rscript $BIN/rexpfit.r ${inputChr} ${jOutput} ${output_col} ${lval} ${hval} ${affine} ${plot} "TRUE" 
	done
	
	snpname=`$BIN/grabpars -p $parfile -x "snpname:" `
	perl $BIN/wtjack.pl -i ${input} -m ${snpname} -b ${BIN}
fi



