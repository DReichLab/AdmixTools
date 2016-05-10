DOCUMENTATION OF F4 ratio estimation (qpF4ratio):

F4 ratio estimation allows inference of the mixing proportions of an admixture event, even without access to accurate surrogates for the ancestral populations.

For any 5 populations that are related to each other as shown in Figure 4 (Ancient Admixture, Patterson et al. 2012), then we can compute the admixture proportion \alpha as - 

\alpha = f4(A,O; X,C)/ f4(A,O; B, C) 

qpF4ratio requires that the input data is available in EIGENSTRAT format.  To convert to the appropriate format, one can use CONVERTF program. See README.CONVERTF for documentation of programs for converting file formats.

Executable and source code:
------------------------------------------------------------------------------
 
For information about installing the program, see README.ADMIXTOOLS. After installing the programs, the executable for F4 ratio estimation (qpF4ratio) should be located in the bin directory.

To run qpF4ratio, type the following on a linux machine. 
$DIR/bin/qpF4ratio -p parfile >logfile

$DIR: Path to the bin directory.
logfile: Name of the logfile. The logfile contains the output of the run. 
parfile: Name of parameter file

DESCRIPTION OF EACH PARAMETER in parfile:

genotypename:   input genotype file (in eigenstrat format)
snpname:   input snp file      (in eigenstrat format)
indivname:   input indiv file    (in eigenstrat format)
popfilename:  list (contains list of 5 poulations- A 0 : X C:: A 0 : B C. Colon separators are not required and are ignored by the program. They have been included to make the input/output more readable). 
blgsize:	jackknife block size (in centimorgan)

Note:
If you would like to use SNP number to define block size, you can:
make a fake .snp file with (false) physical positions, and 0 genetic distances.  
The software defaults to 1M bases = 1cm, so if (for example) you set the
distance between each snp as 10000 bases and blgsize: .01 then each block
will be 100 snps.

DESCRIPTION OF OUTPUT FILE:
The program will write all the output to stdout. The output file prints the parfile entered by the user, number of snps and individuals, jackknife block size, number of blocks for jackknif
e and the results.

The results have the following format -
result:   Pop1 (A)  Pop2 (O) Pop3 (X)  Pop4 (C): Pop1 (A)  Pop2 (O) Pop5 (B)  Pop4 (C)  alpha  std. err  Z (null=0)

The result for each set of 5 populations is shown on a separate line.


See example shown in-
examples/F4ratio.log 

Nick Patterson
<nickp@broadinstitute.org>
------------------------------------------------------------------------------



