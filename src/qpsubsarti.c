#include "qpsubs.h" 
#include "eigsubs.h" 
extern int fancynorm, verbose, plotmode, outnum ;
extern FILE *fstdetails ;

void weightjackfourier(double *est,double *sig,double *mean,double *jmean,double *jwt,int g,double* prho); 
static Indiv **indm ;
static void wjackestx(double *est, double *sig, double mean, double *jmean, double *jwt, int g)  ;
static void wjackvestx(double *vest, double *var, int d, double *mean, double **jmean, double *jwt, int g)  ;


void printsc(int tpat[3][4], double tscore[3], char **eglist, double ymin)  
{
 int a, b, c, d ;
 int *tp, k ;
 
 tp = tpat[0] ;
 printf("qscore: ") ;
 a = tp[0] ; printf ("%15s ", eglist[a]) ;
 a = tp[1] ; printf ("%15s ", eglist[a]) ;
 printf("  ") ;
 a = tp[2] ; printf ("%15s ", eglist[a]) ;
 a = tp[3] ; printf ("%15s ", eglist[a]) ;
 for (k=0; k<3; k++) {
  tp = tpat[k] ;
  printf("%2d ", tp[0]) ;
  printf("%2d ", tp[1]) ;
  printf(" ") ;
  printf("%2d ", tp[2]) ;
  printf("%2d ", tp[3]) ;
  printf("%9.3f", tscore[k]) ;
 }

 printf ("  %9.3f\n", ymin) ;
 printnl() ;

}
void xcopy(int rp[4], int a , int b, int c, int d)  
{

  rp[0] = a ; 
  rp[1] = b ;  
  rp[2] = c ;  
  rp[3] = d ;


}
void
settsc(int tpat[3][4], double tscore[3], int rpat[3][4], double rscore[3]) 
/// process rscore and return scores in tscore with tscore[0] best
{
  double ww[3], w2[3], y ;
  double xmax, xmin ;
  int indx[3], i, k ;

  vvt(ww, rscore, rscore, 3) ;  
  vmaxmin(ww, 3, &xmax, &xmin) ;
  vsp(ww, ww, -xmin, 3) ;  
  copyarr(ww, w2, 3) ; 
  vst(ww, ww, -1.0, 3) ;
  sortit(w2, indx, 3) ;  
  y = w2[1] ;  // second best score
  vsp(ww, ww, y, 3) ; 

  for (i=0; i<3; i++)  { 
   k = indx[i] ; 
   if (i==0) y = rscore[k] ;
   tscore[i] = ww[k] ;
   copyiarr(rpat[k], tpat[i], 4) ;
  }
}

void getpdata(int *rawcol, double *pm, double *pn, int *xtypes, int nrows, int numeg)
{
  int *ytypes, n=0 ;
  int i, k, t, g ;
  double *data, y ;
  
  vzero(pm, numeg) ; 
  vzero(pn, numeg) ;

  ZALLOC(ytypes, nrows, int) ; 
  ZALLOC(data, nrows, double) ; 
   for (k=0; k<nrows; ++k) {  
    g = rawcol[k] ;
    t = xtypes[k] ;
    if (g<0) continue ; 
    if (t<0) continue ; 
    data[n] = g ; 
    ytypes[n] = t ; 
    ++n ; 
   }

   if (n<=1) { 
    free(ytypes) ;
    free(data) ;
    return ;
   }
   y = asum(data, n) / (double) n ;  
// vsp(data, data, -y, n) ; 

   y = 0.5*y ;  
   y = y*(1.0-y) ;

   if (y<.001) {  
    free(ytypes) ;
    free(data) ;
    return ;
   }

   vst(data, data, 1.0/sqrt(y), n) ;
    for (i=0; i<n; i++) {  
     t = ytypes[i] ;
     if (t<0) continue ; 
     if (t>=numeg) continue ;
     pm[t] += data[i] ;
     pn[t] += 1.0  ;
    }

    vsp(pn, pn, 1.0e-8, numeg) ;
    vvd(pm, pm, pn, numeg) ;

   free(ytypes) ;
   free(data) ;


}

void
gethscore(double *hscore, double *scores, 
  int a, int b, int c, int d, int numeg) 

{
 hscore[0] = qhdiff(scores,  a, b, c, d, numeg) ;
 hscore[1] = qhdiff(scores,  a, b, c, d, numeg) ;
 hscore[2] = qhdiff(scores,  a, b, c, d, numeg) ;
}

void
getrscore(double *rscore, double *rho, double **zz, 
  int ncols, int a, int b, int c, int d, int numeg,  int *blabels, int nblocks) 

{
 rscore[0] = qcorr(zz, &rho[0], ncols, a, b, c, d, numeg, blabels, nblocks) ;
 rscore[1] = qcorr(zz, &rho[1], ncols, a, c, b, d, numeg, blabels, nblocks) ;
 rscore[2] = qcorr(zz, &rho[2], ncols, a, d, b, c, numeg, blabels, nblocks) ;
}
double qhdiff(double *scores,  int a, int b, int c, int d, int numeg) 
{
  double tt[4], xmax, xmin ; 
  tt[0] = scores[a*numeg+c] ;
  tt[1] = scores[a*numeg+d] ;
  tt[2] = scores[b*numeg+c] ;
  tt[3] = scores[b*numeg+d] ;
  vmaxmin(tt, 4, &xmax, &xmin) ; 
  return -(xmax-xmin) ;
}

double qcorr(double **zz, double *rho, int ncols, int a, int b, int c, int d, int numeg,  int *blabels, int nblocks) 
{
  double *z1, *z2, y, xrho, xsig ;  
  int u, v ;  

  u = MIN(a, b) ;
  v = MAX(a, b) ; 
  z1 = zz[u*numeg+v] ;

  u = MIN(c, d) ;
  v = MAX(c, d) ; 
  z2 = zz[u*numeg+v] ;

  corrwjack(&xrho, &xsig, z1, z2, ncols, blabels, nblocks) ; 
  *rho = xrho ;
  y = xrho/xsig ;

  return y ;

}
int
loadindx(Indiv **xindlist, int *xindex, Indiv **indivmarkers, int numindivs) 
{
  int i, n=0 ;  
  Indiv *indx ;
  for (i=0; i<numindivs; i++)  {
   indx = indivmarkers[i] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus == NO) continue ;
   xindex[n] = i ;
   xindlist[n] = indx ;
   ++n ;
  }
  return n ;
}
int
loadsnpx(SNP **xsnplist, SNP **snpmarkers, int numsnps, Indiv **indivmarkers) 
{
  int i, n=0 ;  
  SNP *cupt ;
  for (i=0; i<numsnps; i++)  {
   cupt = snpmarkers[i] ;
   if (cupt -> ignore) continue ;
   if (numvalidgt(indivmarkers, cupt) == 0) continue ;
   xsnplist[n] = cupt ;
   ++n ;
  }
  return n ;
}

void
getrawcol(int *rawcol, SNP *cupt, int *xindex, int nrows) 
{
  int t, j ;
  for (j=0; j< nrows; j++) {  
   t = xindex[j] ;
   rawcol[j] = getgtypes(cupt, t) ;          
// if (verbose) printf("www %d %d %d\n", j, t, rawcol[j]) ;
  }
}
void
getrawcolx(int **cc, SNP *cupt, int *xindex, int nrows, Indiv **indm) 
{
  int t, j, g, tt ;
  int *gg ;
  Indiv *indx ;
  static int ncall = 0 ;
  
  ++ncall ;
//  tt = strcmp(cupt -> ID, "rs10914979") ;
  tt = -1 ;
  for (j=0; j< nrows; j++) {  
   t = xindex[j] ;
   gg = cc[j] ;
   ivclear(gg, -1, 2) ;
   g = getgtypes(cupt, t) ;
   if (tt==0) printf("zzcolx %d %d %d\n", j, t, g) ;

   if (ncall==-1) {
     printf("zzindx2:  %s\n", indm[230] -> egroup) ;
     printf("zz1 %d %d %d\n",j, t, g) ;
     indx = indm[t] ;
     printf("yy2 %20s %20s %20s %d %d %d\n", cupt ->ID, indx -> ID, indx -> egroup, j, t, g) ;
   }

   if (g<0) continue ;  
   gg[0] = g ; gg[1] = 2-g ;  
   if (cupt -> chrom != 23) continue ;
   if (indm[t] -> gender != 'M') continue ;  
   if (g==1) { 
    ivclear(gg, -1, 2) ;
    continue ;
   }
   g = g/2 ; 
   gg[0] = g ; gg[1] = 1-g ;  
  }
}


void
getcolx(double *xcol, SNP *cupt, int *xindex, int nrows, int col,  
 double *xmean, double *xfancy)             
// side effect set xmean xfancy
{
 Indiv *indx  ;
 int  j,  n, g, t ;
 double y, pmean, p ;
 int *rawcol ;

  ZALLOC(rawcol, nrows, int) ;
  n = cupt -> ngtypes ;
  if (n<nrows) fatalx("bad snp: %s %d\n", cupt -> ID, n) ;
  getrawcol(rawcol, cupt, xindex, nrows) ;
  floatit(xcol, rawcol, nrows) ;

  vadjust(xcol, nrows, &pmean) ;
  if (fancynorm) {  
   p = 0.5*pmean ;   // autosomes
   y = p * (1.0-p) ;  
   if (y<=0.0) return ;
   y = 1.0/sqrt(y) ;
   vst(xcol, xcol, y, nrows) ;
  }
  else y = 1.0 ;
  if (xmean != NULL) {
   xmean[col] = pmean*y ; 
   xfancy[col] = y ;
  }
  free(rawcol) ;
}
void
loadxdataind(double *xrow, SNP **snplist, int ind,  int ncols)             
{
 SNP *cupt ;
 Indiv *indx ;
 int i, j, k, n,  g ;

 for (i=0; i<ncols; i++) {  
   cupt = snplist[i] ;
   g = getgtypes(cupt, ind) ; 
   xrow[i] = (double) g ;
 }
}
void fixxrow(double *xrow, double *xmean, double *xfancy, int len) 
{
    int i ;

    vvt(xrow, xrow, xfancy, len) ;
    for (i=0; i<len; i++) { 
     if (xrow[i] < -0.1) xrow[i] = 0.0 ; 
     else xrow[i] -= xmean[i] ;
    }
}

void dofancy(double *cc, int n, double *fancy) 
{
  int i, t, nmiss=0 ;
  int top, bot ;  
  double p, yvar, y ;

  top = bot = 0 ;
  for (i=0; i<n; i++) {  
   t = nnint(cc[i]) ;
   if (t<0) { 
    ++nmiss ;
    continue ;
   }
   top += t ;
   bot += 2 ;
  }
  if (bot==0) return ;
  if (top == 0) return ;
  if (top == bot) return ;
  p = (double) top / ((double) bot) ;
  yvar = p*(1.0-p) ;
  y = 1.0/sqrt(yvar) ;
  vst(cc, cc, y, n )  ;                    
  *fancy = y ;
}

int vadjust(double *cc, int n, double *pmean) 
/* take off mean  force missing to zero */
{
 double ynum, ysum, y, ymean ;
 int i, nmiss=0 ;

 ynum = ysum = 0.0 ;
 for (i=0; i<n; i++) {  
  y = cc[i] ;
  if (y < 0.0) { 
   ++nmiss ;
   continue ;
  }
  ++ynum ;
  ysum += y ;
 }
 if (ynum==0.0) fatalx("logic bug all missing\n") ;
 ymean = ysum/ynum ;
 for (i=0; i<n; i++) {  
  y = cc[i] ;
  if (y < 0.0) cc[i] = 0.0 ;
  else cc[i] -= ymean ;
 }
 if (pmean != NULL) *pmean = ymean ;
 return nmiss ;
}

double yll(double x1, double x2, double xlen) 
{
  double m1, m2, var ;

  if (xlen < 1.5 ) return 0.0 ;
  m1 = x1/xlen ;
  m2 = x2/xlen ;  
  var = m2-m1*m1 ;
  if (var <= 0.0) fatalx("bad yll\n") ;
  return -0.5*xlen*log(var) ;
}
void
calcmean(double *wmean, double *vec, int len, int *xtypes, int numeg) 
{

   int i, k ;
   double y1 ;
   double *w0, *popsize ;

   ZALLOC(w0, len, double) ;
   ZALLOC(popsize, numeg, double) ;

   y1 = asum(vec, len)/ (double) len ;  // mean
   vsp(w0, vec, -y1, len) ;

    for (i=0; i<len; i++)  { 
     k = xtypes[i] ;
     ++popsize[k] ;
     wmean[k] += w0[i] ;
    }



    vsp(popsize, popsize, 1.0e-12, numeg) ;
    vvd(wmean, wmean, popsize, numeg) ;

    free(w0) ;
    free(popsize) ;



}

void
setmiss(SNP **snpm, int numsnps) 
{
   SNP *cupt ;
   int i, j, t, n, tot ;

   for (i=0; i<numsnps; i++)  {  
    cupt = snpm[i] ;
    n = cupt -> ngtypes ;
    if (n <= 0) continue ;
    tot = 0 ;
    for (j=0; j<n; j++) {  
     if (getgtypes (cupt, j) >= 0) { 
      t = 1 ;
     }
     else { 
      t = 0 ;
     }
     putgtypes(cupt, j, t) ;
     tot +=t ;
    }
    if (verbose)
     printf("Valids: %s %d\n", cupt -> ID, tot) ;
   }
}

void
setfvecs(double *fvecs, double *evecs, int nrows, int numeigs) 
// plotmode each eigenvector min 0 max 1
{

   double *w ;
   double xmax, xmin ;
   int i, j ;

   ZALLOC(w, nrows, double) ;

   for (j=0; j<numeigs; j++) {  
    copyarr(evecs+j*nrows, w, nrows) ;
    if (plotmode==NO) {   
     vst(fvecs+j*nrows, w, 10.0, nrows) ;
     continue ;
    }
    copyarr(w, fvecs+j*nrows, nrows) ;
   }
   free(w) ;
}
void countpops(int ***counts, SNP **xsnplist, int *xindex, int *xtypes, int nrows, int ncols) 
// countpops is int [ncols][npops][2]  already zero filled
{
     int col, i, g1, g2, g, k1 ; 
     SNP *cupt ;
     int *rawcol ;

     ZALLOC(rawcol, nrows, int) ;
     for (col = 0; col < ncols; ++col) {  
      cupt = xsnplist[col] ;
      getrawcol(rawcol, cupt, xindex, nrows) ;
      for (i=0; i<nrows; i++) { 
       g = rawcol[i] ;
       k1 = xtypes[i] ;
       if (k1<0) continue ;
       g1 = 0 ; 
       if (k1>0) g1 = 1 ;  
       g2 = g-g1 ;
       ++counts[col][k1][g1] ;
       ++counts[col][k1][g2] ;
      }
     }
     free(rawcol) ;
}

// setidmat used to scale

void fixrho (double *a, int n) 
// turn a into correlation matrix
{
   double *d, *tt, y ;
   int i ;  

   ZALLOC(d, n, double) ;
   ZALLOC(tt, n*n, double) ;

   getdiag(d, a, n) ;

   vsqrt(d, d, n) ;
   addouter(tt, d, n) ;

   vvd(a, a, tt, n*n) ;


  free(d) ;
  free(tt) ;


}
void printdiag(double *a, int n) 
{
   double *d, *tt, y ;
   int i ;  

   ZALLOC(d, n, double) ;
   getdiag(d, a, n) ;
    y = asum(d,n) / (double) n ; 
    for (i=0; i<n; i++) {  
     printf("diag: %9.3f\n", d[i]/y) ;
    }


   free(d) ;
   abort() ;

}
int
ridoutlier(double *evecs, int n, int neigs, 
 double thresh, int *badlist, OUTLINFO **outinfo) 
{
/* badlist contains list of outliers */
 double *ww, y1 , y2 ; 
 int *vbad ;
 int i, j  ;
 int nbad = 0 ; 
 OUTLINFO *outpt;

 ZALLOC(ww, n, double) ;
 ZALLOC(vbad, n, int) ;
 for(j=0;j<n;j++)  {
   outpt = outinfo[j];
   outpt->vecno = -1;
 }
 for (i=0; i<neigs; ++i) {  
  copyarr(evecs+i*n, ww, n) ;
  y1 = asum(ww, n) / (double) n ;
  vsp(ww, ww, -y1, n) ;
  y2 = asum2(ww, n) / (double)  n ;
  y2 = sqrt(y2) ;
  vst(ww, ww, 1.0/y2, n) ;

  for (j=0; j<n; j++) {  
   if (fabs(ww[j])>thresh) { 
    vbad[j] = 1 ;
    outpt = outinfo[j];
    if (outpt->vecno < 0)  {
      outpt->vecno = i;
      outpt->score = ww[j];
    }
   }
  }
 }
 for (j=0; j<n; j++) {  
  if (vbad[j] == 1) { 
   badlist[nbad] = j ;
   ++nbad ;
  }
 }
 free(ww) ; 
 free(vbad) ;
 return nbad ;

}

void
addoutersym(double *X, double *v, int n)  
{
  int i, j ;

  for (i=0; i<n; i++) { 
    for (j=i; j<n; j++) {
      X[i*n+j] += v[i]*v[j] ;
    }
  }
}

void symit(double *X, int n) 
{
  int i, j ;

  for (i=0; i<n; i++) { 
    for (j=i+1; j<n; j++) {
      X[j*n+i] = X[i*n+j]  ;
    }
  }
}

double divcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) 
/* heterozygosity for 2 pops */
{
   int c1[2], c2[2], *cc ;
   int *rawcol ;
   int k, g, i ; 
   double ya, yb, yaa, ybb, p1, p2, en, ed ;
   double z, zz, h1, h2, yt ;


   ZALLOC(rawcol, nrows, int) ;

   getrawcol(rawcol, cupt, xindex, nrows)  ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;
    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (cc == NULL) continue ;
    g = rawcol[i] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    cc[1] += 2-g ;
   }

   ya = c1[0] ;
   yb = c1[1] ;
   yaa = c2[0] ;
   ybb = c2[1] ;
   z = ya + yb ;
   zz = yaa+ybb ;
   if ((z<0.1) || (zz<0.1)) { 
    *estn = 0.0 ;
    *estd = -1.0; /* no data */
    free(rawcol) ;
    return 0.0;
   }

   en = ya*ybb + yb*yaa ; 
   ed = z*zz ;

   *estn = en ; 
   *estd = ed ; 
   

   free(rawcol) ;
   return z + zz ;

}

void f3y(double *estn,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3) 
{
   int c1[2], c2[2], c3[2], *cc ;
   int *rawcol ;
   int k, g, i, a, b  ; 
   double ya, yb, yaa, ybb, p1, p2, p3, en, ed ;
   double z, zz, h1, h2, yt ;
   double ywt ; 


   ZALLOC(rawcol, nrows, int) ;

   getrawcol(rawcol, cupt, xindex, nrows)  ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;
   ivzero(c3, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;
    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (k==type3) cc = c3 ; 
    if (cc == NULL) continue ;
    g = rawcol[i] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    cc[1] += 2-g ;
   }

   ya = a = c1[0] ;
   yb = b = c1[1] ;
   z = ya + yb ;


   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   yaa = c2[0] ;
   ybb = c2[1] ;
   yt = yaa+ybb ;
   p2 = yaa/yt ;         

   yaa = c3[0] ;
   ybb = c3[1] ;
   yt = yaa+ybb ;
   p3 = yaa/yt ;         

   en = (p1-p2)*(p1-p3) ;  
   en -= h1/z ; 
 
   *estn = en ; 
   

   free(rawcol) ;

}

void f2sc(double *estn,  double *estd, SNP *cupt, Indiv **indm,  
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3) 
// processes X chromosome correctly
{
   int c1[2], c2[2], c3[2], c4[2], *cc ;
   int *rawcol ;
   int k, g, i, a, b  ; 
   double ya, yb, yaa, ybb, p1, p2, p3, p4, en, ed ;
   double z, zz, h1, h2, h3, yt ;
   double z2, z3 ;
   double ywt ; 
   int **ccc, *ccpt[3] ;


   ccc = initarray_2Dint(nrows, 2, 0) ;


   getrawcolx(ccc, cupt, xindex, nrows, indm) ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;
   ivzero(c3, 2) ;

   ccpt[0] = c1 ;
   ccpt[1] = c2 ;
   ccpt[2] = c3 ;

   *estn = 0 ;  
   *estd = -1 ;

   for (i=0; i< nrows; i++)  { 

    k = xtypes[i] ;
    cc = NULL ;

    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (k==type3) cc = c3 ; 

    if (cc == NULL) continue ;

    g = ccc[i][0] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    g = ccc[i][1] ;
    cc[1] += g ;

   }


/**
   printf("qq1: %d  ", cupt -> markernum) ;
   printimat(c1, 1, 2) ;
*/

   for (i=0; i<=2; i++) { 
    cc = ccpt[i] ;
    a = intsum(cc,2) ;
    if (a<2) {  
     free2Dint(&ccc, nrows) ; 
     return ;
    }
   }

   ya = a = c1[0] ;
   yb = b = c1[1] ;
   z = ya + yb ;


   yt = ya+yb ;
   p1 = ya/yt ;       

   h1 = ya*yb/(yt*(yt-1.0)) ;

   yaa = c2[0] ;
   ybb = c2[1] ;
   z2 = yt = yaa+ybb ;
   h2 = yaa*ybb/(yt*(yt-1.0)) ;
   p2 = yaa/yt ;         

   yaa = c3[0] ;
   ybb = c3[1] ;
   z3 = yt = yaa+ybb ;
   h3 = yaa*ybb/(yt*(yt-1.0)) ;
   p3 = yaa/yt ;         

// h1 0 is OK trap if necessary in calling program

   en = (p2-p3)*(p2-p3) ;
   en -= h2/z2 ;           
   en -= h3/z3 ;           

   if (isnan(en)) fatalx("f3 bug\n") ;
 
   *estn = en ; 
   *estd = 2.0*h1 ;
   

   free2Dint(&ccc, nrows) ;

}
void f3sc(double *estn,  double *estd, SNP *cupt, Indiv **indm,  
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3) 
// processes X chromosome correctly
{
   int c1[2], c2[2], c3[2], c4[2], *cc ;
   int *rawcol ;
   int k, g, i, a, b  ; 
   double ya, yb, yaa, ybb, p1, p2, p3, p4, en, ed ;
   double z, zz, h1,  yt ;
   double ywt ; 
   int **ccc, *ccpt[3] ;


   ccc = initarray_2Dint(nrows, 2, 0) ;


   getrawcolx(ccc, cupt, xindex, nrows, indm) ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;
   ivzero(c3, 2) ;

   ccpt[0] = c1 ;
   ccpt[1] = c2 ;
   ccpt[2] = c3 ;

   *estn = 0 ;  
   *estd = -1 ;

   for (i=0; i< nrows; i++)  { 

    k = xtypes[i] ;
    cc = NULL ;

    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (k==type3) cc = c3 ; 

    if (cc == NULL) continue ;

    g = ccc[i][0] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    g = ccc[i][1] ;
    cc[1] += g ;

   }


/**
   printf("qq1: %d  ", cupt -> markernum) ;
   printimat(c1, 1, 2) ;
*/

   for (i=0; i<=2; i++) { 
    cc = ccpt[i] ;
    a = intsum(cc,2) ;
    if (a<2) {  
     free2Dint(&ccc, nrows) ; 
     return ;
    }
   }

   ya = a = c1[0] ;
   yb = b = c1[1] ;
   z = ya + yb ;


   yt = ya+yb ;
   p1 = ya/yt ;       

   h1 = ya*yb/(yt*(yt-1.0)) ;

   yaa = c2[0] ;
   ybb = c2[1] ;
   yt = yaa+ybb ;
   p2 = yaa/yt ;         

   yaa = c3[0] ;
   ybb = c3[1] ;
   yt = yaa+ybb ;
   p3 = yaa/yt ;         

// h1 0 is OK trap if necessary in calling program

   en = (p1-p2)*(p1-p3) ;
   en -= h1/z ;           

   if (isnan(en)) fatalx("f3 bug\n") ;
 
   *estn = en ; 
   *estd = 2.0*h1 ;
   

   free2Dint(&ccc, nrows) ;

}

void f4yx(double *estn,  SNP *cupt, Indiv **indm,  
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3, int type4) 
// processes X chromosome correctly
{
   int c1[2], c2[2], c3[2], c4[2], *cc ;
   int *rawcol ;
   int k, g, i, a, b  ; 
   double ya, yb, yaa, ybb, p1, p2, p3, p4, en, ed ;
   double z, zz, h1, h2, yt ;
   double ywt ; 
   int **ccc ;


   ccc = initarray_2Dint(nrows, 2, 0) ;


   getrawcolx(ccc, cupt, xindex, nrows, indm) ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;
   ivzero(c3, 2) ;
   ivzero(c4, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;

    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (k==type3) cc = c3 ; 
    if (k==type4) cc = c4 ; 

    if (cc == NULL) continue ;

    g = ccc[i][0] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    g = ccc[i][1] ;
    cc[1] += g ;
   }

   ya = a = c1[0] ;
   yb = b = c1[1] ;
   z = ya + yb ;


   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   yaa = c2[0] ;
   ybb = c2[1] ;
   yt = yaa+ybb ;
   p2 = yaa/yt ;         

   yaa = c3[0] ;
   ybb = c3[1] ;
   yt = yaa+ybb ;
   p3 = yaa/yt ;         

   yaa = c4[0] ;
   ybb = c4[1] ;
   yt = yaa+ybb ;
   p4 = yaa/yt ;         
   en = (p1-p2)*(p3-p4) ;  
 
   *estn = en ; 
   
   free2Dint(&ccc, nrows) ;

}


void f4y(double *estn,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3, int type4) 
{
   int c1[2], c2[2], c3[2], c4[2], *cc ;
   int *rawcol ;
   int k, g, i, a, b  ; 
   double ya, yb, yaa, ybb, p1, p2, p3, p4, en, ed ;
   double z, zz, h1, h2, yt ;
   double ywt ; 


   ZALLOC(rawcol, nrows, int) ;

   getrawcol(rawcol, cupt, xindex, nrows)  ;

   ivzero(c1, 2) ;
   ivzero(c2, 2) ;
   ivzero(c3, 2) ;
   ivzero(c4, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;
    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (k==type3) cc = c3 ; 
    if (k==type4) cc = c4 ; 
    if (cc == NULL) continue ;
    g = rawcol[i] ;
    if (g<0) continue ;  
    cc[0] += g ; 
    cc[1] += 2-g ;
   }

   ya = a = c1[0] ;
   yb = b = c1[1] ;
   z = ya + yb ;


   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   yaa = c2[0] ;
   ybb = c2[1] ;
   yt = yaa+ybb ;
   p2 = yaa/yt ;         

   yaa = c3[0] ;
   ybb = c3[1] ;
   yt = yaa+ybb ;
   p3 = yaa/yt ;         

   yaa = c4[0] ;
   ybb = c4[1] ;
   yt = yaa+ybb ;
   p4 = yaa/yt ;         
   en = (p1-p2)*(p3-p4) ;  
 
   *estn = en ; 
   

   free(rawcol) ;

}



void fstcolyy(double *estnmat, double *estdmat, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int numeg) 
/**
  NP style n, d estimation for fst No ascertainment  
 like fstcoly but a matrix of populations so data is only accessed once 
*/
{
   int *c1, *c2, *cc ;
   int *rawcol ;
   int k, g, i, j, a, b  ; 
   double ya, yb, yaa, ybb, p1, p2, en, ed ;
   double z, zz, h1, h2, yt ;
   double ywt ; 
   int **ccc, *gg, **ddd ;
   static int ncall = 0 ;


   ++ncall ;
   ccc = initarray_2Dint(nrows, 2, 0) ;
   ddd = initarray_2Dint(numeg, 2, 0) ;

   

   vzero(estnmat, numeg*numeg) ;
   vclear(estdmat, -1.0, numeg*numeg) ;

   for (a=0; a<numeg; a++) { 
    estdmat[a*numeg+a] = 0.0 ;
   }

   if (indm == NULL) {
    ZALLOC(rawcol, nrows, int) ;
    getrawcol(rawcol, cupt, xindex, nrows)  ;
    for (a=0; a<nrows; a++)  {
     g = rawcol[a] ;  
     ccc[a][0] = g ;
     ccc[a][1] = 2-g ;
    }
    free(rawcol) ;
   }

   else {
    getrawcolx(ccc, cupt, xindex, nrows, indm)  ;
   }


   ywt = 1.0 ;  

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;

    if (k<0) continue ;
    if (k>=numeg) continue ;

    cc = ddd[k] ;
    gg = ccc[i] ;
    g = gg[0] ;
    if (g<0) continue ;  
    ivvp(cc, cc, gg, 2) ;
   }

   for (i=0; i<numeg; i++) { 
    for (j=i+1; j<numeg; j++) { 
     c1 = ddd[i] ;
     c2 = ddd[j] ;
     ya = a = c1[0] ;
     yb = b = c1[1] ;
     yaa = c2[0] ;
     ybb = c2[1] ;
     zz = yaa+ybb ;
     z = ya + yb ;
     if ((z<1.5) || (zz<1.5)) {
      continue ;
     }


      z = ya+yb ;

      yt = ya+yb ;
      p1 = ya/yt ;       
      h1 = ya*yb/(yt*(yt-1.0)) ;  // 2 h1 is heterozygosity

      yt = yaa+ybb ;
      p2 = yaa/yt ;         
      h2 = yaa*ybb/(yt*(yt-1.0)) ;

      en = (p1-p2)*(p1-p2) ;  
      en -= h1/z ; 
      en -= h2/zz ; 
 
      ed = en ; 
      ed += h1 ; 
      ed += h2 ;

      if (ed<0.0) fatalx("logic bug\n") ;
      estnmat[i*numeg+j] = estnmat[j*numeg+i] = en*ywt ;
      estdmat[i*numeg+j] = estdmat[j*numeg+i] = ed*ywt ;
    }
   }

   free2Dint(&ccc, nrows) ;
   free2Dint(&ddd, numeg) ;

}



double fstcoly(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) 
/** NP style n, d estimation for fst No ascertainment  */
{
   int c1[2], c2[2], *cc ;
   int *rawcol ;
   int k, g, i, a, b  ; 
   double ya, yb, yaa, ybb, p1, p2, en, ed ;
   double z, zz, h1, h2, yt ;
   double ywt ; 
   int **ccc, *gg ;
   static int ncall = 0 ;


   ++ncall ;
   ccc = initarray_2Dint(nrows, 2, 0) ;

   if (indm == NULL) {
    ZALLOC(rawcol, nrows, int) ;
    getrawcol(rawcol, cupt, xindex, nrows)  ;
    for (a=0; a<nrows; a++)  {
     g = rawcol[a] ;
     ccc[a][0] = g ;
     ccc[a][1] = 2-g ;
    }
    free(rawcol) ;
   }

   else {
    getrawcolx(ccc, cupt, xindex, nrows, indm)  ;
   }


   ivzero(c1, 2) ;
   ivzero(c2, 2) ;

   for (i=0; i< nrows; i++)  { 
    k = xtypes[i] ;
    cc = NULL ;
    gg = ccc[i] ;
    if (ncall==-11) {       
     printf("zzindx1:  %s\n", indm[230] -> egroup) ;
     printf("zz2 %d %d ", type1, type2) ;
     printf("%3d %d %3d %3d\n", i, k, gg[0], gg[1]) ;
    }
    if (k==type1) cc = c1 ; 
    if (k==type2) cc = c2 ; 
    if (cc == NULL) continue ;
    g = gg[0] ;
    if (g<0) continue ;  
    ivvp(cc, cc, gg, 2) ;
   }

   ya = a = c1[0] ;
   yb = b = c1[1] ;
   yaa = c2[0] ;
   ybb = c2[1] ;
   zz = yaa+ybb ;
   z = ya + yb ;
   if ((z<1.5) || (zz<1.5)) {
    *estn =  0.0 ;
    *estd = -1.0 ; /* no data in column */
    free2Dint(&ccc, nrows) ;
    return 0.0;
   }

   ywt = ya*yb/(z*(z-1.0)) ; // z must be at least 2 
   ywt = 1.0 ;  

   z = ya+yb ;

   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;  // 2 h1 is heterozygosity

   yt = yaa+ybb ;
   p2 = yaa/yt ;         
   h2 = yaa*ybb/(yt*(yt-1.0)) ;

   en = (p1-p2)*(p1-p2) ;  
   en -= h1/z ; 
   en -= h2/zz ; 
 
   ed = en ; 
   ed += h1 ; 
   ed += h2 ;

   if (ed<0.0) fatalx("logic bug\n") ;

   *estn = en*ywt ; 
   *estd = ed*ywt ; 
  
/**
   printf("zz %20s %2d %2d  ", cupt ->ID, type1, type2) ;
   printf("%3d %3d ", c1[0], c1[1]) ;
   printf("%3d %3d ", c2[0], c2[1]) ;

   printf(" %9.3f %9.3f", *estn, *estd) ;
   printnl() ;
*/

   free2Dint(&ccc, nrows) ;
   return z + zz ;

}
void
 setplimit(Indiv **indivmarkers, int numindivs, 
 char **eglist, int numeg, int plimit)  
{
     int *indnums ;
     int *psize ;
     int i, k, kk ;
     Indiv *indx ; 

     ZALLOC(indnums, numindivs, int) ;
     ZALLOC(psize, numeg, int) ;


     idperm(indnums, numindivs) ;
     ranperm(indnums, numindivs) ;

     for (i=0; i<numindivs; i++)  {  
      k = indnums[i] ; 
      indx = indivmarkers[k] ;
      if (indx -> ignore) continue ;
      kk = indxindex(eglist, numeg, indx -> egroup) ;
      if (kk<0) continue ; 
      ++psize[kk] ; 
      if (psize[kk] > plimit) indx -> ignore = YES ;
     }



     free(psize) ;
     free(indnums) ;

}

double dohzg(double *top, double *bot, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg)  

{

   int t1, t2 ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int k, g, i, col, j ; 
   double ya, yb, y ;
   double *xtop, *xbot ;
   SNP *cupt ;
   

   vzero(top, numeg*numeg) ;
   vzero(bot, numeg*numeg) ;

   ZALLOC(rawcol, nrows, int) ;
   ZALLOC(pop0, numeg, int) ;
   ZALLOC(pop1, numeg, int) ;
   ZALLOC(popall, numeg, int) ;

   for (col=0; col<ncols;  ++col)  {
    ivzero(popall, numeg) ;
    ivzero(pop0, numeg) ;
    ivzero(pop1, numeg) ;
    cupt = xsnplist[col] ;
    getrawcol(rawcol, cupt, xindex, nrows)  ;
    for (i=0; i< nrows; i++)  { 
     k = xtypes[i] ;
     g = rawcol[i] ;
     if (g<0) continue ;  
     pop1[k] += g ; 
     pop0[k]  += 2-g ;
     popall[k] += 2 ;  // code needs chamging for X  
    }
    for (k=0; k<numeg; k++) {  
     ya = pop0[k] ;
     yb = pop1[k] ;
     top[k*numeg+k] += 2*ya*yb ;
     y = ya + yb ;
     bot[k*numeg+k] += y*(y-1.0) ;
     for (j=k+1; j<numeg; j++) {  
      ya = pop0[j] ;
      yb = pop1[k] ;
      y = ya + yb ;
      top[k*numeg+j] += ya*yb ;
      ya = pop1[j] ;
      yb = pop0[k] ;
      top[j*numeg+k] = top[k*numeg+j] += ya*yb ;

      ya = popall[k] ; 
      yb = popall[j] ; 
      bot[k*numeg+j] += ya*yb ;      

      top[j*numeg+k] = top[k*numeg+j] ;         
      bot[j*numeg+k] = bot[k*numeg+j] ;         
     }
    }
  }
  ZALLOC(xtop, numeg*numeg, double) ;
  ZALLOC(xbot, numeg*numeg, double) ;
  copyarr(bot, xbot, numeg*numeg) ;
  y = bal1(xbot, numeg*numeg) ;
  vst(xtop, top, 1.0/y, numeg*numeg)   ;

  free(xtop) ;
  free(xbot) ;


    free(rawcol) ; 
    free(pop0) ; 
    free(pop1) ; 
    free(popall) ;

}
void setblocks(int *block, int *bsize, int *nblock, SNP **snpm, int numsnps, double blocklen)  
// block, bsize are first element and block length 
// must have been allocated
{
  int n = 0, i ;  
  int chrom, xsize, lchrom, olds ;  
  double fpos, dis, gpos ;  
  SNP *cupt ;
      

  lchrom = -1 ; xsize = 0 ;

  fpos = -1.0e20 ;
  for (i=0; i<numsnps; i++) {  
   cupt = snpm[i] ;
   cupt -> tagnumber = -1 ;
   if (cupt -> ignore) continue ; 
   if (cupt -> isfake) continue ;
   chrom = cupt -> chrom ; 
   gpos = cupt -> genpos ;
   dis = gpos - fpos ;
   if ((chrom != lchrom) || (dis>=blocklen)) {
    if (xsize>0) {  
     block[n] = olds  ;  
     bsize[n] = xsize ;  
     ++n ;  
    }
     lchrom = chrom ; 
     fpos = gpos ; 
     olds = i ;
     xsize = 0 ;
   }
   cupt -> tagnumber = n ;
   ++xsize ;
  }
    if (xsize>0) {  
     block[n] = olds  ;  
     bsize[n] = xsize ;  
     ++n ;  
    }
    *nblock = n ;
    

    return ;
}   

int numblocks(SNP **snpm, int numsnps, double blocklen)  
{
  int n = 0, i ;  
  int chrom, xsize, lchrom, olds ;  
  double fpos, dis, gpos ;  
  SNP *cupt ;
      

  lchrom = -1 ; xsize = 0; olds = -1 ;

  fpos = 0 ;
  for (i=0; i<numsnps; i++) {  
// printf("zz %d %d %d %d\n", i, n, xsize, olds) ;
   cupt = snpm[i] ;
   cupt -> tagnumber = -1 ;
   if (cupt -> ignore) continue ; 
   if (cupt -> isfake) continue ;
   chrom = cupt -> chrom ; 
   gpos = cupt -> genpos ;
   dis = gpos - fpos ;
   if ((chrom != lchrom) || (dis>blocklen)) {
    if (xsize>0) {  
     cupt = snpm[olds] ; 
     printf("info abt block: %d  %d %d\n", n, cupt->chrom, (int) cupt -> physpos) ;
     ++n ;  
    }
     lchrom = chrom ; 
     fpos = gpos ; 
     olds = i ;
     xsize = 1 ;
     continue ;
   }
   ++xsize ;
  }
  return n+1 ;
}   

void corrwjack(double *xrho, double *xsig, double *z1, double *z2, int ncols, int *bcols, int nblocks)
{
   double *gdot, *dot, *wdot ; 
   double **bdot ;
   double *djack, *wjack ;
   double rho, jest, jsig ;
   double y1, y2 ;
   int bnum, i, k ;


   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   ZALLOC(gdot, 6, double) ;
   ZALLOC(wdot, 6, double) ;

   bdot  = initarray_2Ddouble(nblocks, 6,  0.0) ;


   for (i=0; i<ncols; i++) {
    bnum = bcols[i] ;          
    if (bnum<0) continue ;  
    ++wjack[bnum] ;
    dot = bdot[bnum] ;  
    y1 = z1[i] ; 
    y2 = z2[i] ;  
    dot[0] += y1*y1 ;
    dot[1] += y2*y2 ;
    dot[2] += y1*y2 ;
    dot[3] += y1 ; 
    dot[4] += y2 ; 
    dot[5] += 1.0 ;
   }
   for (k=0; k<nblocks; k++) {  
    dot = bdot[k] ;
    vvp(gdot, gdot, dot, 6) ;
   }
   rho = crho(gdot) ;
// printmatw(gdot, 1, 6, 6) ;
   for (k=0; k<nblocks; k++) {  
    dot = bdot[k] ;
    vvm(wdot, gdot, dot, 6) ;
    djack[k] = crho(wdot) ;
   }
   wjackest(&jest, &jsig, rho, djack, wjack, nblocks) ;
   *xrho = jest ;
   *xsig = jsig ;
   
   free(djack) ; 
   free(wjack) ; 
   free(gdot) ; 
   free(wdot) ; 
   free2D(&bdot, nblocks) ;


}
double crho(double *stats) 
{
/* correlation from 6 sufficient statistics */
 double m1, m2,  top, bot, b1, b2, rr ;
 double s1, s2, s11, s22, s12, yn ;
 static int ncall = 0 ;

 ++ncall ;
 s11 = stats[0] ;
 s22 = stats[1] ;
 s12 = stats[2] ;
 s1 = stats[3] ;
 s2 = stats[4] ;
 yn = stats[5] ; 

 m1 = s1/yn ; 
 m2 = s2/yn ;
 top = s12 - yn*m1*m2 ;
 b1 = s11 - yn*m1*m1 ;
 b2 = s22 - yn*m2*m2 ;

 if (ncall < -1) { 
  printf("%9.3f\n", m1) ;
  printf("%9.3f\n", m2) ;
  printf("%9.3f\n", top) ;
  printf("%9.3f\n", b1) ;
  printf("%9.3f\n", b2) ;
 }
 rr = top/sqrt(b1*b2) ; 
 
 return rr ;
}

void setbcols(SNP **xsnplist, int ncols, int *bcols) 
{
   int col, bnum ;  
   SNP *cupt ;
   
   ivclear(bcols, -1, ncols) ;
   for (col=0; col<ncols;  ++col)  {
    cupt = xsnplist[col] ;
    bnum = cupt -> tagnumber ; 
    bcols[col] = bnum ;
   }
}

double
doadmlin(double *jest, double *jsig, double *zlin, double *var, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale, Indiv **indm)     
{

   int t1, t2, kret ;
   int a, b, c ;
   int ng3, ng2 ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int k, g, i, col, j, d ; 
   double ya, yb, y, mean ;
   SNP *cupt ;
   double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop ;
   double **btop, **bbot, wt ;
   double *w1, *w2, *w3 ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   double xest, xsig, ynominal ;
   int bnum  ;

   double *f3, *f3sig  ;
   double *estmat, *zl, *rhs, errest ;
   double *vmean, **vjmean ;
   
   ng2 = numeg*numeg ;
   ng3 = numeg*numeg*numeg ;

   ZALLOC(f3, ng3, double) ;
   ZALLOC(f3sig, ng3, double) ;
   ZALLOC(w1, ng3+2, double) ;
   ZALLOC(w2, ng3+2, double) ;
   ZALLOC(estmat, ng3, double) ;
   ZALLOC(w3, ng3, double) ;
   ZALLOC(gtop, ng3, double) ;
   ZALLOC(gbot, ng3, double) ;
   ZALLOC(wtop, ng3, double) ;
   ZALLOC(wbot, ng3, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   btop  = initarray_2Ddouble(nblocks, ng3,  0.0) ;
   bbot  = initarray_2Ddouble(nblocks, ng3,  0.0) ;

   d = numeg - 1 ;
   vjmean  = initarray_2Ddouble(nblocks, numeg,  0.0) ;
   ZALLOC(vmean, numeg, double) ;

   zl = w1 ; 
   rhs = w2 ;  // overloading

   for (col=0; col<ncols;  ++col)  {
    cupt = xsnplist[col] ;
    if (cupt -> ignore) continue ;
    wt = cupt -> weight ;  
    if (wt <= 0.0) continue ;
    bnum = cupt -> tagnumber ; 
    if (bnum<0) continue ;
    ++wjack[bnum] ;
    top = btop[bnum] ; 
    bot = bbot[bnum] ;

    kret = f3yyx(estmat,  cupt, xindex, xtypes, nrows, numeg, indm) ;
    if (kret < 0) continue ;
    vst(estmat, estmat, wt*scale, ng3) ;
    vvp(top, top, estmat, ng3) ;
    vsp(bot, bot, 1.0, ng3) ;

   }

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, ng3) ;
     vvp(gbot, gbot, bot, ng3) ;
    }

   vsp(w2, gbot, 1.0e-10, ng3) ;
   vvd(f3, gtop, w2, ng3) ;

   vzero(zl, numeg) ;
   estmix(zl+1, f3, numeg) ;
   copyarr(zl+1, vmean, d) ;


  ynominal = y = estmix(zlin, f3, numeg) ;

  if (verbose) {

   for (i=0; i<numeg; ++i) { 
    printf("f3: base number %d:\n", i) ;
    printmatw(f3+i*numeg*numeg, numeg, numeg, numeg) ;
   }

   printf("nominal error: %12.6f\n", y) ;
  }



   ytop = ybot = errest = 0.0 ;

   vvd(wtop, gtop, gbot, ng3) ;  // delete-block estimate

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvm(wtop, gtop, top, ng3) ;
     vvm(wbot, gbot, bot, ng3) ;
     vsp(wbot, wbot, 1.0e-10, ng3) ;
     vvd(wtop, wtop, wbot, ng3) ;  // delete-block estimate
     vzero(zl, numeg) ;
     djack[k] = estmix(zl+1, wtop, numeg) ;
     copyarr(zl+1, vjmean[k], d) ;
///  printf("yyy: %4d  %9.3f  %12.6f\n", k, wjack[k], djack[k]) ;
     mulmat(rhs, top, zl, numeg, numeg, 1) ;
     y = vdot(zl, rhs, numeg) ;
     ytop += y ;
     ybot += bot[0] ;
     if (verbose) 
      printf("www: %4d  %9.3f  %12.6f\n", k, wjack[k], y) ;
    }

    errest = ytop/ybot ;
// jackknife estimate of standard error for variance
    wjackest(&xest, &xsig, ynominal, djack, wjack, nblocks) ;
    wjackvest(vmean, var,  d, zlin, vjmean, wjack, nblocks)  ;
    *jest = xest ;
    *jsig = xsig ;

    free(w1) ; 
    free(w2) ; 
    free(w3) ; 
    free(estmat) ;

    free(gbot) ; 
    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;
    free(f3) ;
    free(f3sig) ;

    free(vmean) ;

    free2D(&btop, nblocks);
    free2D(&bbot, nblocks);
    free2D(&vjmean, nblocks);

    return errest ;

}


void
dof3(double *f3, double *f3sig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale, int mode)     
{

   int t1, t2 ;
   int a, b, c ;
   int ng3 ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   SNP *cupt ;
   double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop ;
   double **btop, **bbot, wt ;
   double *w1, *w2, *w3 ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   int bnum  ;
   double *estmat ;
   
   ng3 = numeg*numeg*numeg ;
   ZALLOC(w1, ng3, double) ;
   ZALLOC(w2, ng3, double) ;
   ZALLOC(estmat, ng3, double) ;
   ZALLOC(w3, ng3, double) ;
   ZALLOC(gtop, ng3, double) ;
   ZALLOC(gbot, ng3, double) ;
   ZALLOC(wtop, ng3, double) ;
   ZALLOC(wbot, ng3, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   btop  = initarray_2Ddouble(nblocks, ng3,  0.0) ;
   bbot  = initarray_2Ddouble(nblocks, ng3,  0.0) ;

   for (col=0; col<ncols;  ++col)  {
    cupt = xsnplist[col] ;
    if (cupt -> ignore) continue ;
    wt = cupt -> weight ;  
    if (wt <= 0.0) continue ;
    bnum = cupt -> tagnumber ; 
    if (bnum<0) continue ;
    ++wjack[bnum] ;
    top = btop[bnum] ; 
    bot = bbot[bnum] ;

    f3yy(estmat,  cupt, xindex, xtypes, nrows, numeg) ;

    if (mode != 2) {
     vst(estmat, estmat, wt, ng3) ;
     vvp(top, top, estmat, ng3) ;
     vsp(bot, bot, 1.0, ng3) ;
    }
    else {  
     vvp(top, top, estmat, ng3) ;
     vsp(bot, bot, 1.0/wt, ng3) ;
    }
   }

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, ng3) ;
     vvp(gbot, gbot, bot, ng3) ;
    }

   vsp(w2, gbot, 1.0e-10, ng3) ;
   vvd(f3, gtop, w2, ng3) ;


    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvm(wtop, gtop, top, ng3) ;
     vvm(wbot, gbot, bot, ng3) ;
     vsp(wbot, wbot, 1.0e-10, ng3) ;
     vvd(top, wtop, wbot, ng3) ;  // delete-block estimate
    }
    vsp(gbot, gbot, 1.0e-10, ng3) ;
    vvd(gtop, gtop, gbot, ng3) ;
    
    
    for (a=0; a<numeg ; a++) { 
     for (b=0; b<numeg ; b++) { 
      for (c=0 ; c<numeg ; c++) { 
       if (a==b) continue ;
       if (a==c) continue ;
       if (c<b)  continue ;
      for (k=0; k<nblocks; k++) {  
       top = btop[k] ; 
       djack[k] = dump3(top, a, b, c, numeg) ;
      }
      
      mean = dump3(gtop, a, b, c, numeg) ;
      wjackest(&jest, &jsig, mean, djack, wjack, nblocks) ;
      bump3(f3sig, a, b, c, numeg, jsig) ;
      bump3(f3sig, a, c, b, numeg, jsig) ;
     }
    }
   }
   vst(f3, f3, scale, ng3) ;
   vst(f3sig, f3sig, scale, ng3) ;

    free(w1) ; 
    free(w2) ; 
    free(w3) ; 
    free(estmat) ;

    free(gbot) ; 
    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;

    free2D(&btop, nblocks);
    free2D(&bbot, nblocks);

}
void bump2(double *x, int a, int b, int n, double val)  
{
   int k ;
   k = a ; 
   k *= n ; 
   k += b ; 
   x[k] += val ;
}
void bump3(double *x, int a, int b, int c, int n, double val)  
{
   int k ;
   k = a ; 
   k *= n ; 
   k += b ; 
   k *= n ;
   k += c ; 
   x[k] += val ;
}
double dump2(double *x, int a, int b, int n)  
{
   int k ;
   double val ;
   k = a ; 
   k *= n ; 
   k += b ; 
   val = x[k] ;
   return val ;
}
double dump3(double *x, int a, int b, int c, int n)  
{
   int k ;
   double val ;
   k = a ; 
   k *= n ; 
   k += b ; 
   k *= n ;
   k += c ; 
   val = x[k] ;
   return val ;
}
void
dof4(double *f4, double *f4sig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale, int mode)     
{

   int t1, t2 ;
   int a, b, c, d ;
   int ng4 ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   SNP *cupt ;
   double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop ;
   double **btop, **bbot, wt ;
   double *w1, *w2, *w3 ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   int bnum  ;
   int nloop = 0 ;
   
   ng4 = numeg*numeg*numeg*numeg ;
   ZALLOC(w1, ng4, double) ;
   ZALLOC(w2, ng4, double) ;
   ZALLOC(w3, ng4, double) ;
   ZALLOC(gtop, ng4, double) ;
   ZALLOC(gbot, ng4, double) ;
   ZALLOC(wtop, ng4, double) ;
   ZALLOC(wbot, ng4, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   btop  = initarray_2Ddouble(nblocks, ng4,  0.0) ;
   bbot  = initarray_2Ddouble(nblocks, ng4,  0.0) ;

   for (col=0; col<ncols;  ++col)  {
    cupt = xsnplist[col] ;
    if (cupt -> ignore) continue ;
    wt = cupt -> weight ;  
    if (wt <= 0.0) continue ;
    bnum = cupt -> tagnumber ; 
    if (bnum<0) continue ;
    if (bnum>=nblocks) fatalx("logic bug\n") ;
    ++wjack[bnum] ;
    top = btop[bnum] ; 
    bot = bbot[bnum] ;

    for (a=0; a<numeg ; a++) { 
     for (b=0; b<numeg ; b++)  { 
      for (c=0 ; c<numeg ; c++)  { 
       for (d=0 ; d<numeg ; d++)  { 

       if (a==b) continue ;
       if (a==c) continue ;
       if (a==d) continue ;
       if (b==c) continue ;
       if (b==d) continue ;
       if (c==d) continue ;

       if (b<a)  continue ;
       if (c<a)  continue ;
       if (d<a)  continue ;
       if (d<c)  continue ;

       f4y(&ytop,  cupt, xindex, xtypes, nrows, a, b, c, d) ;
       ++nloop ;  
      //  if (nloop<100) printf("zz1 %d %d %d %d %9.3f\n", a, b, c, d, ytop)  ;
       if (isnan(ytop)) fatalx("zznan\n") ;

       if (mode != 2) {
        bump4x(top, a, b, c, d, numeg, wt*ytop) ;
        bump4x(top, b, a, c, d, numeg, -wt*ytop) ;
        bump4x(bot, a, b, c, d, numeg, 1.0) ;
        bump4x(bot, b, a, c, d, numeg, 1.0) ;
       }
       else {
        bump4x(top, a, b, c, d, numeg, ytop) ;
        bump4x(top, b, a, c, d, numeg, -ytop) ;
        bump4x(bot, a, b, c, d, numeg, 1.0/wt) ;
        bump4x(bot, b, a, c, d, numeg, 1.0/wt) ;
       }
       
      }
     }
    }
   }
  }   

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, ng4) ;
     vvp(gbot, gbot, bot, ng4) ;
    }

   vsp(w2, gbot, 1.0e-10, ng4) ;
   vvd(f4, gtop, w2, ng4) ;


    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvm(wtop, gtop, top, ng4) ;
     vvm(wbot, gbot, bot, ng4) ;
     vsp(wbot, wbot, 1.0e-10, ng4) ;
     vvd(top, wtop, wbot, ng4) ;  // delete-block estimate
    }
    vsp(gbot, gbot, 1.0e-10, ng4) ;
    vvd(gtop, gtop, gbot, ng4) ;
    
    
    for (a=0; a<numeg ; a++) { 
     for (b=0; b<numeg ; b++) { 
      for (c=0 ; c<numeg ; c++) { 
       for (d=0 ; d<numeg ; d++) { 
        if (a==b) continue ;
        if (a==c) continue ;
        if (a==d) continue ;
        if (b==c) continue ;
        if (b==d) continue ;
        if (c==d) continue ;

        if (b<a)  continue ;
        if (c<a)  continue ;
        if (d<a)  continue ;
        if (d<c)  continue ;

       for (k=0; k<nblocks; k++) {  
        top = btop[k] ; 
        djack[k] = dump4(top, a, b, c, d, numeg) ;
       }
      
       mean = dump4(gtop, a, b, c, d, numeg) ;
       wjackest(&jest, &jsig, mean, djack, wjack, nblocks) ;
       bump4x(f4sig, a, b, c, d, numeg, jsig) ;
       bump4x(f4sig, b, a, c, d, numeg, jsig) ;
      }
     }
    }
   }
   vst(f4, f4, scale, ng4) ;
   vst(f4sig, f4sig, scale, ng4) ;

    free(w1) ; 
    free(w2) ; 
    free(w3) ; 

    free(gbot) ; 
    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;

    free2D(&btop, nblocks);
    free2D(&bbot, nblocks);

}

void bump4x(double *x, int a, int b, int c, int d, int n, double val)  
{
   bump4(x, a, b, c, d, n, val) ;
   bump4(x, b, a, d, c, n, val) ;
   bump4(x, c, d, a, b, n, val) ;
   bump4(x, d, c, b, a, n, val) ;
}

void bump4(double *x, int a, int b, int c, int d, int n, double val)  
{
   int k ;
   k = a ; 
   k *= n ; 
   k += b ; 
   k *= n ;
   k += c ; 
   k *= n ; 
   k += d ;
   x[k] += val ;
}
void set4x(double *x, int a, int b, int c, int d, int n, double val)  
{
   set4(x, a, b, c, d, n, val) ;
   set4(x, b, a, d, c, n, val) ;
   set4(x, c, d, a, b, n, val) ;
   set4(x, d, c, b, a, n, val) ;
}

void set4(double *x, int a, int b, int c, int d, int n, double val)  
{
   int k ;
   k = a ; 
   k *= n ; 
   k += b ; 
   k *= n ;
   k += c ; 
   k *= n ; 
   k += d ;
   x[k] = val ;
}

double dump4(double *x, int a, int b, int c, int d, int n)  
{
   int k ;
   double val ;
   k = a ; 
   k *= n ; 
   k += b ; 
   k *= n ;
   k += c ; 
   k *= n ; 
   k += d ;
   val = x[k] ;
   return val ;
}

double
dofstnumx(double *fst, double *fstest, double *fstsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, Indiv **indivmarkers, int fstmode)     

// fstmode is classic mode (smartpca)
// fstmode 2  is fstdmode

{

   int t1, t2 ;
   int a, b ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int t, k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   SNP *cupt ;
   double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop ;
   double **btop, **bbot, wt, rho ;
   double *w1, *w2, *w3 ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   int bnum  ;
   int nloop = 0, fstdnum =0  ;
   double *ztop, *zbot, qtop, qbot ;
   char **eglist ;  
   int jfourierflag = 1 ;

   indm = indivmarkers ;
   if ((fstdetails != NULL) && (indm == NULL)) fatalx("bug in dofstnumx\n") ;

   ZALLOC(eglist, numeg, char *) ;
   for (k=0; k< nrows; ++k) {  
    if (indm == NULL) break ;
    j = xtypes[k] ; 
    if (j<0) continue ; 
    if (j>=numeg) continue ; 
    t = xindex[k] ; 
    eglist[j] = indm[t] -> egroup ;
   }
   
   ZALLOC(w1, numeg*numeg, double) ;
   ZALLOC(w2, numeg*numeg, double) ;
   ZALLOC(w3, numeg*numeg, double) ;
   ZALLOC(gtop, numeg*numeg, double) ;
   ZALLOC(gbot, numeg*numeg, double) ;
   ZALLOC(wtop, numeg*numeg, double) ;
   ZALLOC(wbot, numeg*numeg, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   ZALLOC(ztop, numeg*numeg, double) ;
   ZALLOC(zbot, numeg*numeg, double) ;
   btop  = initarray_2Ddouble(nblocks, numeg*numeg,  0.0) ;
   bbot  = initarray_2Ddouble(nblocks, numeg*numeg,  0.0) ;

   vzero(fst, numeg*numeg) ;
   vzero(fstest, numeg*numeg) ;
   vzero(fstsig, numeg*numeg) ;


   for (col=0; col<ncols;  ++col)  {
    cupt = xsnplist[col] ;
    if (cupt -> ignore) continue ;
    wt = cupt -> weight ;  
    if (wt <= 0.0) continue ;
    bnum = cupt -> tagnumber ; 
    if (bnum<0) continue ;
    ++wjack[bnum] ;
    top = btop[bnum] ; 
    bot = bbot[bnum] ;

    fstcolyy(ztop, zbot, cupt, xindex, xtypes, nrows, numeg) ;

    for (a=0; a<numeg ; a++) { 
     for (b=a+1; b<numeg ; b++) { 
      k = a*numeg+b ; 
      ytop = ztop[k] ;
      ybot = zbot[k] ;
     if (fstdetails != NULL) {  
      if (fstdnum==0) { 
       fprintf(fstdetails,"%15s ", "## pop 1") ;
       fprintf(fstdetails,"%15s ", "pop 2") ;
       fprintf(fstdetails,"%15s ", "snpname") ;  
       fprintf(fstdetails,"%12s ", "N") ;
       fprintf(fstdetails,"%12s ", "D") ;
       fprintf(fstdetails, "\n") ;
      }
      fprintf(fstdetails,"%15s ", eglist[a]) ;
      fprintf(fstdetails,"%15s ", eglist[b]) ;
      fprintf(fstdetails,"%15s ", cupt -> ID) ; 
      fprintf(fstdetails,"%12.6f ", ytop) ;
      fprintf(fstdetails,"%12.6f ", ybot) ;
      fprintf(fstdetails, "\n") ;
      ++fstdnum ;
     }


      if (ybot<0.0) continue ;


      if (fstmode == NO) {
       top[k] += wt*ytop ; 
       bot[k] += 1.0  ;
      }

      if (fstmode == YES) {
       top[k] += ytop ; 
       bot[k] += ybot ;
      }

      if (fstmode == 2) {
       top[k] += ytop ; 
       bot[k] += 1.0/wt ;
      }

      w1[k] += ytop ;
      w2[k] += ybot ;
// classic fst estimate

     }
    }
   }
// symmetrize
   for (a=0; a<numeg ; a++) { 
    for (b=a+1; b<numeg ; b++) { 
      top[b*numeg+a] = top[a*numeg+b] ;
      bot[b*numeg+a] = bot[a*numeg+b] ;
      w1[b*numeg+a]  = w1[a*numeg+b] ;
      w2[b*numeg+a]  = w2[a*numeg+b] ;
    }
   }
   
// printf("zzz ") ; printmat(wjack, 1, nblocks) ;

   vsp(w2, w2, 1.0e-10, numeg*numeg) ;
   vvd(fst, w1, w2, numeg*numeg) ;
   

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, numeg*numeg) ;
     vvp(gbot, gbot, bot, numeg*numeg) ;
    }

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvm(wtop, gtop, top, numeg*numeg) ;
     vvm(wbot, gbot, bot, numeg*numeg) ;
     vsp(wbot, wbot, 1.0e-10, numeg*numeg) ;
     vvd(top, wtop, wbot, numeg*numeg) ;  // delete-block estimate
    }
    vsp(gbot, gbot, 1.0e-10, numeg*numeg) ;
    vvd(gtop, gtop, gbot, numeg*numeg) ;
    
    
    for (i=0; i<numeg; i++)  {  
     for (j=i+1; j<numeg ; j++) {  
      for (k=0; k<nblocks; k++) {  
       top = btop[k] ; 
       djack[k] = top[i*numeg+j] ;
      }
      
      ++nloop ;
      mean = gtop[i*numeg+j] ;
      
      printf("Mean, nblocks: %f %d\n", mean,nblocks);

      if(jfourierflag == 1)
	{
	  weightjackfourier(&jest, &jsig, &mean, djack, wjack, nblocks,&rho) ;
	  printf("rho, mean from fourier jk:%2.6f %f\n",rho,mean);
	}
      else
	{
	  wjackest(&jest, &jsig, mean, djack, wjack, nblocks) ;
	}

      fstest[i*numeg+j] = fstest[j*numeg+i] = jest ;
      fstsig[i*numeg+j] = fstsig[j*numeg+i] = jsig ;

       if (nloop == -1)  {
         printf("ddd\n") ;
         printf("mean: %9.3f\n", mean) ;
         printmat(djack, 1, nblocks) ;
         printmat(wjack, 1, nblocks) ;
         printf("%9.3f %9.3f\n", jest, jsig) ;
       }
     }
    }

/**
   printf("fst:\n") ;
   printmat(fst, numeg, numeg) ;
   printnl() ;
   printmat(fstest, numeg, numeg) ;
   printnl() ;
   printmat(fstsig, numeg, numeg) ;
   printnl() ;
*/

    yscal = 1.0 ;
    if (fstmode != YES) {
     copyarr(fstsig, w3, numeg*numeg) ;
     vsp(w3, w3, 1.0e-10, numeg*numeg) ;
     vvd(w1, fst, w3, numeg*numeg) ;
     vvd(w2, fstest, w3, numeg*numeg) ;
//   now do regression  w1 = yscal * w2
     y1 = vdot(w1, w2, numeg*numeg) ;
     y2 = vdot(w2, w2, numeg*numeg) ;
     yscal = y1/y2 ;
     vst(fstest, fstest, yscal, numeg*numeg) ;
     vst(fstsig, fstsig, yscal, numeg*numeg) ;
    }


    free(eglist) ;
    free(w1) ; 
    free(w2) ; 
    free(w3) ; 

    free(gbot) ; 
    free(wtop) ; 
    free(wbot) ; 
    free(ztop) ; 
    free(zbot) ; 
    free(djack) ;
    free(wjack) ;

    free2D(&btop, nblocks);
    free2D(&bbot, nblocks);
    printf("Done with wjackest2\n");

    return yscal ;

}

double
dofstnum(double *fst, double *fstest, double *fstsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks)     
{

   dofstnumx(fst, fstest, fstsig, xsnplist, xindex, xtypes, nrows, ncols, numeg, nblocks, NULL, NO) ;

}

void setmgpos(SNP **snpm, int numsnps, double *maxgdis) 
// find max genetic distance
{
  double minpos, maxdis ; 
  int chrom, lchrom, i ;
  SNP *cupt ;

  minpos = 99999.0 ; 
  lchrom = -1 ;

  maxdis = -9999 ;
  for (i=0; i<numsnps; i++) {  
   cupt = snpm[i] ;
   chrom = cupt -> chrom ;  
   if (chrom != lchrom) { 
    lchrom = chrom ; 
    minpos = cupt -> genpos ;
   }
   maxdis = MAX(maxdis, cupt -> genpos - minpos) ;
  }
  *maxgdis = maxdis ;
}

void setgfromp(SNP **snpm, int numsnps)  
{
  int  i ;
  SNP *cupt ;


  for (i=0; i<numsnps; i++) {  
   cupt = snpm[i] ;
   cupt -> genpos = (cupt -> physpos) / 1.0e8 ;
  }
}

void wjackest(double *est, double *sig, double mean, double *jmean, double *jwt, int g)  
// test for jwt 0 
{
  double *jjmean, *jjwt ;
  int i, n ;

  ZALLOC(jjmean, g, double) ;
  ZALLOC(jjwt, g, double) ;
  n = 0 ;

  for (i=0; i<g ; ++i)  {  
   if (jwt[i] < 1.0e-6) continue ;
   jjmean[n] = jmean[i] ;
   jjwt[n] = jwt[i] ;
   ++n ;
  }

  wjackestx(est, sig, mean, jjmean, jjwt, n) ; 
  free(jjmean) ;
  free(jjwt) ;
}

static void wjackestx(double *est, double *sig, double mean, double *jmean, double *jwt, int g)  
// weighted jackknife see wjack.tex
// mean is natural estimate.  jmean[k] mean with block k removed.  jwt is weight for block (sample size)
{

  double *tdiff, *hh, *xtau, *w1, *w2 ;  
  double jackest, yn, yvar ;  
  int k ;

  if (g<=1) fatalx("(wjackest) number of blocks <= 1\n") ;
  ZALLOC(tdiff, g, double) ;
  ZALLOC(hh, g, double) ;
  ZALLOC(xtau, g, double) ;
  ZALLOC(w1, g, double) ;
  ZALLOC(w2, g, double) ;

  yn = asum(jwt, g) ;  
 
  vsp(tdiff, jmean, -mean, g) ; 
  vst(tdiff, tdiff, -1.0, g) ;  
  jackest = asum(tdiff, g) + vdot(jwt, jmean, g)/yn ;
// this is equation 2

  vclear(hh, yn, g) ;  
  vvd(hh, hh, jwt, g) ;
/**
  for (k=0; k<g; ++k) {
   if (jwt[k] > 0.0) hh[k] /= jwt[k] ;  
   else hh[k] *= 1.0e20 ;
  }
*/
// jwt should be positive

  vst(xtau, hh, mean, g) ; 
  vsp(w1, hh, -1.0, g) ; 
  vvt(w2, w1, jmean, g) ; 
  vvm(xtau, xtau, w2, g) ;
   
  vsp(xtau, xtau, -jackest, g) ;
  vvt (xtau, xtau, xtau, g) ;     
  vvd (xtau, xtau, w1, g) ; 
  yvar = asum(xtau, g) / (double) g ;   
  *est = jackest ;  
  *sig = sqrt(yvar) ;

   free(tdiff) ;
   free(hh) ;
   free(xtau) ;
   free(w1) ;
   free(w2) ;

}

void ndfst5(double *zzest, double *zzsig, double **zn, double **zd, int ncols, int *bcols, int nblocks)
{
#define NPAR  5
   double *djack, *wjack ;
   double qest, jest, jsig ;
   double y1, y2 ;
   int bnum, i, k ;
   int a, b, c ;
   double *gn, *gd, **xn, **xd, *xx, *qqest, *test, *tn, *td, **xqest ;

   ZALLOC(gn, 4*4, double)  ;
   ZALLOC(gd, 4*4, double)  ;
   ZALLOC(tn, 4*4, double)  ;
   ZALLOC(td, 4*4, double)  ;
   ZALLOC(qqest, NPAR, double)  ;

   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;

   xn  = initarray_2Ddouble(nblocks, 4*4,  0.0) ;
   xd  = initarray_2Ddouble(nblocks, 4*4,  0.0) ;
   xqest  = initarray_2Ddouble(nblocks, NPAR,  0.0) ;


   for (i=0; i<ncols; i++) {
    bnum = bcols[i] ;          
    if (bnum<0) continue ;  
    if (bnum>=nblocks) fatalx("bad bug\n") ;
    ++wjack[bnum] ;
    for (a=0; a<4; a++) {  
     for (b=a+1; b<4; b++) {  
      c = 4*a+b ;
      xn[bnum][c] += zn[i][c] ;
      xd[bnum][c] += zd[i][c] ;
     }
    }
   }
   for (k=0; k<nblocks; k++) {  
    xx = xn[k] ;
    vvp(gn, gn, xx, 4*4) ;
    xx = xd[k] ;
    vvp(gd, gd, xx, 4*4) ;
   }
   verbose = YES ;
   regestit(qqest, gn, gd) ;
   printf("qqest: ") ;  
   printmatw(qqest, 1, 5, 5) ;
   verbose = NO ;

   for (k=0; k<nblocks; k++) {  
    xx = xn[k] ;
    vvm(tn, gn, xx, 4*4) ;
    xx = xd[k] ;
    vvm(td, gd, xx, 4*4) ;
    regestit(xqest[k], tn, td) ;
   }
   for (a=0; a<NPAR; a++)   {
    for (k=0; k<nblocks; ++k) { 
     djack[k] = xqest[k][a]  ;
    }
    wjackest(&jest, &jsig, qqest[a], djack, wjack, nblocks) ;
    zzest[a] = jest ; 
    zzsig[a] = jsig ;
   }
   
   free2D(&xqest, nblocks) ;
   free2D(&xn, nblocks) ;
   free2D(&xd, nblocks) ;
   free(djack) ; 
   free(wjack) ; 
   free(gn) ;
   free(gd) ;
   free(tn) ;
   free(td) ;
   free(qqest) ;
}

void regestit(double *ans, double *xn, double *xd)
{
 int a, b, c, k ;  
 double *co, *rr ;
 double f ;

 ZALLOC(co, 6*5, double) ;
 ZALLOC(rr, 6, double) ;

/**
 printf("zzreg\n") ;
 printmat(xn, 4, 4) ;
 printnl() ;
 printmat(xd, 4, 4) ;
 printnl() ;
*/

 verbose = NO ;

 k=0; a=0; b=1; 
 c = 4*a+b ; f = xn[c]/xd[c] ;
 co[k*5+0] = co[k*5+1] = 1 ;  
 rr[k] = f ;

 k=1; a=2; b=3; 
 c = 4*a+b ; f = xn[c]/xd[c] ;
 co[k*5+3] = co[k*5+4] = 1 ;  
 rr[k] = f ;

 k=2; a=0; b=2; 
 c = 4*a+b ; f = xn[c]/xd[c] ;
 co[k*5+0] = co[k*5+2] =  co[k*5+3] = 1 ;  
 rr[k] = f ;

 k=3; a=0; b=3; 
 c = 4*a+b ; f = xn[c]/xd[c] ;
 co[k*5+0] = co[k*5+2] =  co[k*5+4] = 1 ;  
 rr[k] = f ;

 k=4; a=1; b=2; 
 c = 4*a+b ; f = xn[c]/xd[c] ;
 co[k*5+1] = co[k*5+2] =  co[k*5+3] = 1 ;  
 rr[k] = f ;

 k=5; a=1; b=3; 
 c = 4*a+b ; f = xn[c]/xd[c] ;
 co[k*5+1] = co[k*5+2] =  co[k*5+4] = 1 ;  
 rr[k] = f ;

 regressit(ans, co, rr, 6, 5) ;

 free(co) ;
 free(rr) ;

}

void
setwt(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, int nrows, 
  int *xindex, int *xtypes, char *outpop, char **eglist, int numeg) 
{
 int *rawcol ; 
 SNP *cupt ; 
 int i, k, j, t, kk, maxeg ;
 int a0, a1, aa ;
 int **ccx, **ccc, *cc ;
 double wt, p ;
 int a, g ;

 t = strcmp(outpop, "NONE") ;  
 if (t==0) outnum = -1  ;
 maxeg = MAX(outnum, numeg) + 1 ;
 ccx  = initarray_2Dint(maxeg, 2,  0) ;
 ccc  = initarray_2Dint(nrows, 2, 0) ;
 t = -1 ;

// printf("zzqq %d %d\n", outnum, numeg) ;
 for (i=0; i<numsnps; ++i) {  
  cupt = snpmarkers[i] ; 
  cupt -> weight = 0 ;  
//  t = strcmp(cupt -> ID, "rs10914979") ;
  if (cupt -> ignore) continue ;

  getrawcolx(ccc, cupt, xindex, nrows, indivmarkers)  ;
  iclear2D(&ccx, maxeg, 2, 0) ;
  for (k=0; k<nrows; ++k) { 
   a = xtypes[k] ;

   if (i==-1)  {  
    printf("zzq %d %d %d  %d\n", i, k, outnum, ccc[k][0]) ;
   }

   if (a<0) continue ; 
   if (a>=maxeg) continue ;
   g = ccc[k][0] ;
   if (g<0) continue ;
   cc = ccx[a] ;
   ivvp(cc, cc, ccc[k], 2) ;
  }

  if (outnum<0) {  
   a0 = a1 = 0 ;
   for (j=0; j< numeg; ++j) {
    a0 += ccx[j][0] ;
    a1 += ccx[j][1] ;
   }
  }

  else {
   a0 = ccx[outnum][0] ;
   a1 = ccx[outnum][1] ;
  }

  aa = a0 + a1 ;
  if (a0==0) continue ;
  if (a1==0) continue ;
  p = (double) a0 / (double) aa ;  
  wt = 1.0/(p*(1.0-p)) ;
  if (outnum == -99) wt = 1.0 ;

   if (t==0) {
    for (k=0; k<nrows; ++k) {
     printf("ww1: %d %d %d ", k, xtypes[k], xindex[k]) ;
     printimat(ccc[k], 1, 2) ;
    }
   }
  for (k=0; k<numeg ; ++k)  {  
   a0 = ccx[k][0] ;
   a1 = ccx[k][1] ;
   aa = a0+a1 ; 
   if (t==0) printf("zzyy %d %d %d\n", k, a0, a1) ;
   
   if (aa<2) { 
     wt = 0 ;
     break ;
   }
   if (k<numeg) continue ;
  }
  cupt -> weight = wt ;
 }

 for (i=0; i<numsnps; ++i) {  
  cupt = snpmarkers[i] ; 
  if (cupt -> weight <= 0.0) cupt -> ignore = YES ;
 }
 

 free2Dint(&ccx, maxeg) ; 
 free2Dint(&ccc, nrows) ; 

}

void
countg(int *rawcol, int **cc, int *xtypes, int n, int ntypes)  
{
  int g, i, c0, c1, k ;

  iclear2D(&cc, ntypes, 2, 0) ;
  for (i=0; i<n; i++) {  
   g = rawcol[i] ;
   if (g<0) continue ; 
   c0 = g ; 
   c1 = 2-g ; 
   k = xtypes[i] ;  
   if (k<0) continue ; 
   if (k>ntypes) continue ;
   cc[k][0] += c0 ; 
   cc[k][1] += c1 ;  
  }
}   

void
dohzgjack(double *hest, double *hsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int *bcols, int nblocks)     
{

   int t1, t2 ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   SNP *cupt ;
   double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop ;
   double **btop, **bbot ;
   int bnum  ;
   
   ZALLOC(gtop, numeg*numeg, double) ;
   ZALLOC(gbot, numeg*numeg, double) ;
   ZALLOC(wtop, numeg*numeg, double) ;
   ZALLOC(wbot, numeg*numeg, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   btop  = initarray_2Ddouble(nblocks, numeg*numeg,  0.0) ;
   bbot  = initarray_2Ddouble(nblocks, numeg*numeg,  0.0) ;

   ZALLOC(rawcol, nrows, int) ;
   ZALLOC(pop0, numeg, int) ;
   ZALLOC(pop1, numeg, int) ;
   ZALLOC(popall, numeg, int) ;

   ivclear(bcols, -1, ncols) ;

   for (col=0; col<ncols;  ++col)  {
    ivzero(popall, numeg) ;
    ivzero(pop0, numeg) ;
    ivzero(pop1, numeg) ;
    cupt = xsnplist[col] ;
    bnum = cupt -> tagnumber ; 
    bcols[col] = bnum ;
    if (bnum<0) continue ;
    ++wjack[bnum] ;
    top = btop[bnum] ; 
    bot = bbot[bnum] ;
    getrawcol(rawcol, cupt, xindex, nrows)  ;
    for (i=0; i< nrows; i++)  { 
     k = xtypes[i] ;
     if (k<0) continue ; 
     if (k>=numeg) continue ;
     g = rawcol[i] ;
     if (g<0) continue ;  
     pop1[k] += g ; 
     pop0[k]  += 2-g ;
     popall[k] += 2 ;  // code needs chamging for X  
    }
    for (k=0; k<numeg; k++) {  
     ya = pop0[k] ;
     yb = pop1[k] ;
     top[k*numeg+k] += 2*ya*yb ;
     y = ya + yb ;
     bot[k*numeg+k] += y*(y-1.0) ;
     for (j=k+1; j<numeg; j++) {  
      ya = pop0[j] ;
      yb = pop1[k] ;
      y = ya + yb ;
      top[k*numeg+j] += ya*yb ;
      ya = pop1[j] ;
      yb = pop0[k] ;
      top[j*numeg+k] = top[k*numeg+j] += ya*yb ;

      ya = popall[k] ; 
      yb = popall[j] ; 
      bot[k*numeg+j] += ya*yb ;      

      top[j*numeg+k] = top[k*numeg+j] ;         
      bot[j*numeg+k] = bot[k*numeg+j] ;         
     }
    }
  }
    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, numeg*numeg) ;
     vvp(gbot, gbot, bot, numeg*numeg) ;
    }
/**
    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, numeg*numeg) ;
     vvp(gbot, gbot, bot, numeg*numeg) ;
    }
*/
    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvm(wtop, gtop, top, numeg*numeg) ;
     vvm(wbot, gbot, bot, numeg*numeg) ;
     vsp(wbot, wbot, 1.0e-10, numeg*numeg) ;
     vvd(top, wtop, wbot, numeg*numeg) ;  // delete-block estimate
    }
    vsp(gbot, gbot, 1.0e-10, numeg*numeg) ;
    vvd(gtop, gtop, gbot, numeg*numeg) ;
    for (i=0; i<numeg; i++)  {  
     for (j=i; j<numeg ; j++) {  
      for (k=0; k<nblocks; k++) {  
       top = btop[k] ; 
       djack[k] = top[i*numeg+j] ;
      }
      
      mean = gtop[i*numeg+j] ;
      wjackest(&jest, &jsig, mean, djack, wjack, nblocks) ;
//    printf("zz %d %d %12.6f %12.6f\n", i, j, mean, jest) ;
      hest[i*numeg+j] = hest[j*numeg+i] = jest ;
      hsig[i*numeg+j] = hsig[j*numeg+i] = jsig ;
     }
    }

    free(rawcol) ; 
    free(pop0) ; 
    free(pop1) ; 
    free(popall) ;
    free(gtop) ; 
    free(gbot) ; 
    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;

    free2D(&btop, nblocks);
    free2D(&bbot, nblocks);

}

void wjackvest(double *vest, double *var, int d, double *mean, double **jmean, double *jwt, int g)  
// test for jwt 0 
{
  double **jjmean, *jjwt ;
  int i, n ;

  jjmean = initarray_2Ddouble(g, d,  0.0) ;
  ZALLOC(jjwt, g, double) ;

  n = 0 ;

  for (i=0; i<g ; ++i)  {  
   if (jwt[i] < 1.0e-6) continue ;
   copyarr(jmean[i], jjmean[n], d) ;
   jjwt[n] = jwt[i] ;
   ++n ;
  }

  wjackvestx(vest, var, d, mean, jjmean, jjwt, n) ; 

  free2D(&jjmean, g) ;
  free(jjwt) ;
}


static
void wjackvestx(double *vest, double *var, int d, double *mean, double **jmean, double *jwt, int g)  
// weighted jackknife see wjack.tex
// mean is natural estimate.  jmean[k] mean with block k removed.  jwt is weight for block (sample size)
/** 
 mean is d long 
 jjmean is [g][d]  giving jackknifed estimates after deleting each block  
 jwt  is jackknife weight (weight 0 should be OK but is deprecated 

 Output vest is d long (jackknife estimate) 
 var is error variance 
*/
{

  double *xtau, *hh ;  
  double *jackest, yn, yvar ;  
  double *wa ;
  int j, k ;
  double y1, y2 ;

  if (g<=1) fatalx("(wjackvest) not enough blocks\n") ;

  ZALLOC(hh, g, double) ;
  ZALLOC(xtau, d, double) ;
  ZALLOC(wa, d, double) ;

  jackest = vest ;
  
  vzero(var, d*d) ;
  vzero(jackest, d) ;

  yn = asum(jwt, g) ;  
 
  for (k=0; k<g; ++k) {
   vvm(wa, mean, jmean[k], d) ;
   vvp(jackest, jackest, wa, d) ;
   vst(wa, jmean[k], jwt[k]/yn, d) ;
   vvp(jackest, jackest, wa, d) ;
  }
// this is equation 2

  vclear(hh, yn, g) ;  

  for (k=0; k<g; ++k) {

   if (jwt[k] > 0.0) hh[k] /= jwt[k] ;  
   else hh[k] *= 1.0e20 ;

   y1 = hh[k] ;
   vst(xtau, mean, y1, d) ;
   --y1 ;  
   vst(wa, jmean[k], y1, d) ;  
   vvm(xtau, xtau, wa, d) ;
   vvm(xtau, xtau, jackest, d) ;       
   y2 = 1.0/sqrt(y1) ;
   vst(wa, xtau, y2, d) ;
   addouter(var, wa, d) ;
  }
// jwt should be positive

   
  vst(var, var, 1.0 / (double) g, d*d) ;   

   free(hh) ;
   free(xtau) ;
   free(wa) ;

}

int f3yyx(double *estmat,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int numeg, Indiv **indm)
{
   int *c1, *c2, *c3, *cc ;
   int *rawcol ;
   int k, g, i, a, b, c  ; 
   int a0, a1, kret ;
   double ya, yb, yaa, ybb, p1, p2, p3, en, ed ;
   double z, zz, h1, h2, yt ;
   double ywt ; 

   int **ccc, *gg, **ccx ;
   static int ncall = 0 ;


   ++ncall ;
   ccc = initarray_2Dint(nrows, 2, 0) ;
   ccx = initarray_2Dint(numeg+1, 2, 0) ;

   vzero(estmat, numeg*numeg*numeg) ;

    getrawcolx(ccc, cupt, xindex, nrows, indm)  ;

    for (k=0; k<nrows; ++k) { 
     a = xtypes[k] ;
     if (a<0) continue ; 
     if (a>=numeg) continue ;
     g = ccc[k][0] ;
     if (g<0) continue ;
     cc = ccx[a] ;
     ivvp(cc, cc, ccc[k], 2) ;
    }

    kret = 1 ;

    for (a=0; a<numeg ; a++) {
     for (b=0; b<numeg ; b++) {
      for (c=0 ; c<numeg ; c++) {
       if (a==b) continue ;
       if (a==c) continue ;
       if (c<b)  continue ;

       c1 = ccx[a] ;
       c2 = ccx[b] ;
       c3 = ccx[c] ;

   ya =  (double) c1[0] ;
   yb =  (double) c1[1] ;
   z = ya + yb ;


   yt = ya+yb ;
   if (yt<=0)  {
     kret = -1; break ; 
   }
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   

   yaa =  (double) c2[0] ;
   ybb =  (double) c2[1] ;
   yt = yaa+ybb ;
   if (yt<=0)  {
     kret = -1; break ; 
   }
   p2 = yaa/yt ;         
   h2 = yaa*ybb/(yt*(yt-1.0)) ;
   zz = yaa + ybb ;

   yaa =  (double) c3[0] ;
   ybb =  (double) c3[1] ;
   yt = yaa+ybb ;
   if (yt<=0)  {
     kret = -1; break ; 
   }
   p3 = yaa/yt ;         

   en = (p1-p2)*(p1-p3) ;  
   en -= h1/z ; 

   if (b==c) en -= h2/zz ;
 
    bump3(estmat, a, b, c, numeg, en) ;
    if (b!=c) bump3(estmat, a, c, b, numeg, en) ;
   }
  }
 }
   

   free2Dint(&ccc, nrows) ;
   free2Dint(&ccx, numeg+1) ;
   return kret ;

}

void f3yy(double *estmat,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int numeg) 
{
   int *c1, *c2, *c3, *cc ;
   int *rawcol ;
   int k, g, i, a, b, c  ; 
   double ya, yb, yaa, ybb, p1, p2, p3, en, ed ;
   double z, zz, h1, h2, yt ;
   double ywt ; 

   int **ccc, *gg, **ccx ;
   static int ncall = 0 ;


   ++ncall ;
   ccc = initarray_2Dint(nrows, 2, 0) ;
   ccx = initarray_2Dint(numeg, 2, 0) ;

   vzero(estmat, numeg*numeg*numeg) ;

    ZALLOC(rawcol, nrows, int) ;
    getrawcol(rawcol, cupt, xindex, nrows)  ;
    for (a=0; a<nrows; a++)  {
     g = rawcol[a] ;
     ccc[a][0] = g ;
     ccc[a][1] = 2-g ;
    }
    free(rawcol) ;

    for (k=0; k<nrows; ++k) { 
     a = xtypes[k] ;
     if (a<0) continue ; 
     if (a>=numeg) continue ;
     g = ccc[k][0] ;
     if (g<0) continue ;
     cc = ccx[a] ;
     ivvp(cc, cc, ccc[k], 2) ;
    }

    for (a=0; a<numeg ; a++) {
     for (b=0; b<numeg ; b++) {
      for (c=0 ; c<numeg ; c++) {
       if (a==b) continue ;
       if (a==c) continue ;
       if (c<b)  continue ;

       c1 = ccx[a] ;
       c2 = ccx[b] ;
       c3 = ccx[c] ;

   ya =  (double) c1[0] ;
   yb =  (double) c1[1] ;
   z = ya + yb ;


   yt = ya+yb ;
   p1 = ya/yt ;       
   h1 = ya*yb/(yt*(yt-1.0)) ;

   yaa =  (double) c2[0] ;
   ybb =  (double) c2[1] ;
   yt = yaa+ybb ;
   p2 = yaa/yt ;         
   h2 = yaa*ybb/(yt*(yt-1.0)) ;
   zz = yaa + ybb ;

   yaa =  (double) c3[0] ;
   ybb =  (double) c3[1] ;
   yt = yaa+ybb ;
   p3 = yaa/yt ;         

   en = (p1-p2)*(p1-p3) ;  
   en -= h1/z ; 

   if (b==c) en -= h2/zz ;
 
    bump3(estmat, a, b, c, numeg, en) ;
    if (b!=c) bump3(estmat, a, c, b, numeg, en) ;
   }
  }
 }
   

   free2Dint(&ccc, nrows) ;
   free2Dint(&ccx, numeg) ;

}


double estmix(double *z, double *f3, int n) 

/**
 We minimize 
 {f - \sum_i w_i x_i}^2 where f is a target population or individual and x_i are 
 mixing coefficients.  We can write this as in:  
 {\sum_i w_i (f-x_i)}^2 - \lam (\sum_i w_i -1)  
 which yields equations 
 \sum_j X_{ij} w_j = \lam  
 where X are f_3 or f_2 coefficients.
 We solve this by setting lambda = 1 and then normalizing the answer to sum 
 to 1.
*/
{

 int a, b ; 
 int d = n-1 ;
 double *co, *rhs, y, *ww, *w1 ;

 ZALLOC(co, d*d, double) ;
 ZALLOC(ww, d*d, double) ;
 ZALLOC(rhs, d, double) ;
 ZALLOC(w1, d, double) ;

 vclear(rhs, 1.0, d) ;

 for (a=0; a<d; a++) { 
  for (b=a; b<d; b++) { 

   y = dump3(f3, 0, a+1,  b+1, n) ;
// works if a = b as dof3 fixed up
   co[a*d+b] = co[b*d+a] = y  ;

  }
 }
 vclear(w1, 1.0, d) ;
 mulmat(rhs, co, w1, d, d, 1) ;
 mulmat(ww, co, co, d, d, d) ;

 solvit(ww, rhs, d, z) ;

 y = asum(z, d) ;
 if (y==0.0) fatalx("z is zero!\n") ;
 vst(z, z, 1.0/y, d) ;

 mulmat(rhs, co, z, d, d, 1) ;
 y = vdot(z, rhs, d) ;

 free(co) ;
 free(rhs) ;
 free(ww) ;
 free(w1) ;

 return y ;

}
double ff3val(double *ff3, int a, int b, int c, int n)  
{
  double y ;

  y =  dump3(ff3, 0, a, a, n) ; 
  y += dump3(ff3, 0, b, c, n) ; 
  y -= dump3(ff3, 0, b, a, n) ;
  y -= dump3(ff3, 0, c, a, n) ;
  return y ;

}

void weightjackfourier(double *est,double *sig,double *mean,double *jmean,double *jwt,int g,double* prho)
{
  
  double mp,mpr,rhonr,rhodr,rho,mdr,S,cS,jpmean;
  int k,i;
  double *jst, *l,*d,*c,*jp,*wt;

  
  double pi,tmean;

  pi = 2.0*acos(0.0)  ;  

  ZALLOC(jst,g+1,double);
  ZALLOC(l,g,double);
  ZALLOC(d,g,double);
  ZALLOC(c,2*g,double);
  ZALLOC(jp,g,double);
  ZALLOC(wt,g,double);

  for(i=0;i<g;i++)
    {
      wt[i] = jwt[i];
    }

  /**Calculate the mean of the jackknifed means*/
  mp = 0;
  mpr = 0;
  mdr = 0;
  for(k=0;k<g;k++)
    {
      mpr += wt[k]*jmean[k];
      mdr += wt[k];
    }
  mp = mpr/mdr;
  printf("Mp is:%f\n",mp);

  tmean = 0;
  for(k=0;k<g;k++)
    {
      jst[k] = jmean[k] - mp;
      /*      printf("Jst %d %f\n",k,jst[k]);*/
      tmean += jst[k];
    }
  tmean = tmean/g;

  printf("tmean is %f\n",tmean);


  rhonr = 0;
  rhodr = 0;
  rho = 0;
  jst[g] = jst[0];
  for(k=0;k<g;k++)
    {
      rhonr += jst[k]*jst[k+1];
      rhodr += jst[k]*jst[k];
    }


  rho = rhonr/rhodr;
  *prho = rho;

  if(rho < 0)
    {
      printf("Exiting there is an error\n");
      rho = 0;
      return;
    }

  for(k=0;k<g;k++)
    {
      l[k] = 1 + 2*rho*cos((2.0*pi*k)/g);
      d[k] = 1.0/sqrt(l[k]);
    }

  for(i=0;i< 2*g;i++)
    {
      c[i] = 0;
    }

  S = 0;
  for(i=0;i<g;i++)
    {
      c[i] = 0;
      for(k=0;k<g;k++)
	{
	  c[i] += d[k]*cos((2.0*pi*k*i)/g); 
	}
      c[i] = c[i]/g;
      S += c[i];
    }

  /*Using the periodicity of the cosine function*/
  for(i=0;i<g;i++)
    {
      c[g+i] = c[i];
    }

  cS = 1.0/(sqrt(1+2.0*rho));

  jpmean = 0.0;
  for(i=0;i<g;i++)
    {
      jp[i] = 0.0;
      for(k=0;k<g;k++)
	{
	  jp[i] += c[k+i]*jmean[k];
	}
      jp[i] = jp[i]/S;
      jpmean += jp[i];
    }
  
  jpmean = jpmean/g;


  
// overwrite jmean 
// wt not changed -- looks wrong .. 

  for(i = 0; i< g; i++)
    {
      printf("info details: %d %f %f %f\n",i,jmean[i],jst[i],jp[i]);
      jmean[i] = jp[i];
    }

  *mean = jpmean;

  wjackest(est, sig, *mean, jmean, wt, g) ;
  printf("S, cS, jpmean are %f %f %f\n",S,cS,jpmean);
  
  free(jst);
  free(l);
  free(d);
  free(c);
  free(jp);
  free(wt);
	     
}
