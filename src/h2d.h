int mkindh2d (Indiv ** indivmarkers, Indiv *** pindm2, int numindivs) ;
int mkindd2h (Indiv ** indivmarkers, Indiv *** pindm2, int numindivs) ; 

void remaph2d (SNP ** snpmarkers, int numsnps, Indiv ** indivmarkers,
	  Indiv ** indm2, int numindivs, int numind2) ;

void remapd2h (SNP ** snpmarkers, int numsnps, Indiv ** indivmarkers,
	  Indiv ** indm2, int numindivs, int numind2) ;
