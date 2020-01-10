DOCUMENTATION OF qp4diff:

It is often useful to directly test the difference of 2 f4 statistics.  
If we merely run qpDstat and difference 2 results, we cannot correctly compute standard errors. 

qp4diff requires that the input data is available in one of the 5 supported formats. 


Executable and source code:
------------------------------------------------------------------------------
 
For information about installing the program, see README.ADMIXTOOLS. After installing the programs, the executable qp4diff) should be located in the bin directory.

To run qp4diff, type the following on a linux machine. 
$DIR/bin/qp4diff -p parfile >logfile

$DIR: Path to the bin directory.
logfile: Name of the logfile. The logfile contains the output of the run. 
parfile: Name of parameter file

DESCRIPTION OF EACH PARAMETER in parfile:

genotypename:   input genotype file (in eigenstrat format)
snpname:   input snp file      (in eigenstrat format)
indivname:   input indiv file    (in eigenstrat format)
popfilename:  list (contains list of 8 populations- A B : C D :: E F : G H  
Colon separators are not required and are ignored by the program. They have been included to make the input/output more readable). 
Populations can repeat beteeen the left and right quadruplets 
blgsize:	jackknife block size (in centimorgan)

*** 2 optional parameters *** 
allsnps:  YES   
## default.  By default only SNPs with valid f4 for both quads are used.  
firstf4mult:  val 
The first f4 is multiplied by val.  Can be used to test if alpha a_1 = a_2 for a fixed alpha.  

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
result:  A B C D : E F G H  f4diff std. error Z  

Nick Patterson
<nickp@broadinstitute.org>
------------------------------------------------------------------------------



