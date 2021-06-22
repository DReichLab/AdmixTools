qpfstats
--------
         DIR:                 /home/np29/biologyx/v41x/yamdir2/testdir 
          S1:                     qdata1  
         S1X:                     qdata1  
   indivname:                DIR/S1X.ind
     snpname:                 DIR/S1.snp
genotypename:                DIR/S1.geno
 poplistname:                      lista
## must be present.  ne popuation / line.  First population is base
 fstatsoutname:                    fstatsa.txt 
## first line is header.  Must be retained
 allsnps:                          YES 
## default NO  -- dangerous bend
 inbreed:                           NO 
## default (match allsnps:)  
 scale:                             NO 
## default YES -- when fstats are scaled to "match" fst in least squares sense 
*** strongly advuse to specify allsnp: and inbreed: in the paramter file. 

If inbreed: NO the base population (first in poplistname) MUST be a true diploid.   
inbreed: NO when there are pseudo-diploid populations is dubious but may be necessary, 
 especially if there is a population with a single sample.  In this case some f2 statistics 
 can't be evaluated and are set in the output to have huge variance (100).  This means 
 that in qpGraph some tail edges are set 0 -- really are undertermined. 

qpfmv 
-----
 fstatsname:     fstatsa.txt
 popfilename:      f4sslist
## 4 pops / line (as in qpDstat) but f2, f3, f4 can be mixed.  So code A C B C for f3 stat 
     fmvoutname:  fmvq2.out 
 printsd:  YES  ## if you want s.err.  defualt NO. 


New parameters for qpAdm (which should be upwards compatible) 
fstatsname:   <output from qpfstats> 
## if present no need to specify .ind .snp .geno 
numboot:    <# of bootstrap samples + antithetical samples> 
## default 1000 

New parameters for qpWave (which should be upwards compatible) 
fstatsname:   <output from qpfstats> 

New parameters for qpGraph (which should be upwards compatible) 
qpGraph allows to to make qpfstats file 
Example 
1) Make qpfstats  :: ~np29/newadm19/src/jtest1/ww2dir/testdir/qpfs_example/parw2 
D2:          /home/np29/broaddatax/v41  
D1:          D2
S1:           jd1               
indivname:    D1/S1.ind
snpname:      D1/S1.snp 
genotypename: D1/S1.geno
fstatsoutname:  fstatsw2    
allsnps:        YES
## there is an option oldallsnps: YES for compatiblity.  Not recommended. 
loadf3:         YES 
## this is needed to get qpfstats calculated
graphname:       graph1   
diag:           .0001

2) Use qpfstats  
fstatsname:  fstatsw2    
graphname:       graph1   
diag:           .0001

*** qpfstats with allsnps: YES  is a much better way to use all the data than 
the older allsnps mode.   Standard errors are often much reduced and the theory 
is now defensible!

*** there was a bad bug in release 7.0 when qpfstats was run with                 
1) allsnps: YES 
2) inbreed: NO 
3) Samples with pseudo-diploid data.    

In the new release if you do this, f-statistics involving psuedo-populations twice (such aas an f_3 with target pseudo-diploid are  
(correctly) flagged as having very large standard error.  



See ./qpfs.pdf for a brief explanation of the theory. 
