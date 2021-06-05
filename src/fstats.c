double fstat(double *fsindex) 
// loadaa called 
{
   int a, b, c, d ; 
   double p1, p2, p3, p4, yy ; 
     a = fsindex[0] ; 
     b = fsindex[1] ; 
     c = fsindex[2] ; 
     d = fsindex[3] ; 

     if (a==b) return 0 ;
     if (c==d) return 0 ;

    p1 = aafreq[a] ; 
    p2 = aafreq[b] ; 
    p3 = aafreq[c] ; 
    p4 = aafreq[d] ; 

    yy = (p1-p2)*(p3-p4) ;  

    if (a==c)  yy += aaxadd[a] ; 
    if (b==d)  yy += aaxadd[b] ; 
    if (a==d)  yy -= aaxadd[a] ; 
    if (b==c)  yy -= aaxadd[b] ; 

    return yy ; 

}

int
dofstats (double *fsmean, double *fssig, int **fsindex, int nfstats, 
       SNP ** xsnplist, int *xindex, int *xtypes,
       int nrows, int ncols, int numeg, int nblocks, double scale)
{
  double *top, *bot, **btop, **bbot, *wjack , yy, wt ; 
  double *gtop, *gbot, *wmean, *w2, *w3  ; 
  double *jest, *jsig, mean ; 

  int bnum, j, k  ; 
  int ngood = 0 ; 

// pass 1.  Jackknife to get sig

  btop = initarray_2Ddouble(nblocks, nfstats, 0.0) ; 
  bbot = initarray_2Ddouble(nblocks, nfstats, 0.0) ; 
  ZALLOC(wjack, nblocks, double) ; 
  ZALLOC(gtop, nfstats, double) ; 
  ZALLOC(gbot, nfstats, double) ; 

  ZALLOC(wmean, nfstats, double) ; 
  ZALLOC(w2, nfstats, double) ; 
  ZALLOC(w3, nfstats, double) ; 

  ZALLOC(jmean, nfstats, double) ; 
  ZALLOC(jsig, nfstats, double) ; 

  for (col = 0; col < ncols; ++col) {
    cupt = xsnplist[col];
    if (cupt->ignore)
      continue;
    wt = cupt->weight;
    if (wt <= 0.0)
      continue;
    bnum = cupt->tagnumber;
    if (bnum < 0) continue;
    if (bnum>=nblocks) fatalx("logic bug\n") ;

    loadaa (cupt, xindex, xtypes, nrows, numeg);

    ++wjack[bnum];
    ++ngood ; 

    top = btop[bnum];
    bot = bbot[bnum];

    for (k=0; k<nfstats; ++k) { 
     yy = fstat(fsindex[k]) ;
     top[k] += wt*yy ; 
     bot[k] += 1  ; 
    }
  }

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    vvp (gtop, gtop, top, nfstats);
    vvp (gbot, gbot, bot, nfstats);
  }


  vsp (w2, gbot, 1.0e-10, nfstats);
  vvd (wmean, gtop, w2, nfstats);

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    vvm (wtop, gtop, top, nfstats);
    vvm (wbot, gbot, bot, nfstats);
    vsp (wbot, wbot, 1.0e-12, nfstats);
    vvd (top, wtop, wbot,nfstats);	// delete-block estimate
  }

  vsp (gbot, gbot, 1.0e-12, nfstats);
  vvd (gtop, gtop, gbot, nfstats);

// scalar jackknife  
//  weightjack(double *est, double *sig, double mean, double *jmean, double *jwt, int g)  ;
  for (j=0; j<nfstats; ++j) { 
   for (k = 0; k < nblocks; k++) {
    jmean[k] = btop[k][j] ; 
    jwt[k]   = wjack[k] ; 
   } 
   mean = gtop[j] ;  
   weightjack(&jest[j], &jsig[j], mean, jmean, jwt, nblocks) ; 
   printf("jest. pass 1 ") ; 
   printimatx(fsindex[j], 1, 4) ; 
   printf("%12.6f ", mean) ; 
   printf("%12.6f ", jmean[j]) ; 
   printf("%12.6f ", jsig[j]) ; 
   printnl() ; 
  }

  free (wmean);
  free (w2);
  free (w3);

  free (gbot);
  free (wtop);
  free (wbot);
  free (wjack);
  free(jmean, jsig) ;

  free2D (&btop, nblocks);
  free2D (&bbot, nblocks);

  return ngood;

}

