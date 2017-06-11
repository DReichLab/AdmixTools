Hints on using qpGraph
----------------------
1) After writing a new graph use qpreroot -g graph -d junk.dot  to check parsing and that the graph is
what you want.

2) optional modes
a)
inbreed: YES  Use if some populations are pseudo-diploid and all population sizes >= 2  
Without this single samples (like MA1) can bus used, but the edge-length leading to the leaf
is not meaningful.
b)
useallsnps: YES
if coverage is low.  But this is theoretically dubious.
Default is that SNPs are only used if
 not monomorphic and data is present for all populations.
3)
On can fix admixture proportions by specifying
lock  V
where V is a vertex that is admixed.
This is sometimes useful.

4)
For very large graphs the program sometimes has trouble finding the globally best fit.

5) edge lengths of zero suggest that the graph should perhaps be modified.

A -> B of length 0 suggests that maybe B should be ancestral to A

======================================================================
Biggest problem:  How good is the fit?   Arbitrarily good fits can be
obtained by adding enough admixture events.

Nick 6/10/17

