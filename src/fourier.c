void weightjackfourier(double *est,double *sig,double mean,double *kmean,double *jwt,int g,double* prho)
{
  
  double mp,mpr,rhonr,rhodr,rho,mdr,S,cS,jpmean;
  int k,i,l;
  double *jst, *d,*c,*jp,*wt, *wk, *qq;
  double y, y1, y2, ymx, ycx, yg, *jmean, gmean ; 

  
  double pi,tmean;

  pi = 2.0*acos(0.0)  ;  
  yg = (double) g ; 

  ZALLOC(jst,g+1,double);
  ZALLOC(d,g,double);
  ZALLOC(c,2*g,double);
  ZALLOC(jmean,2*g,double);
  ZALLOC(jp,g,double);
  ZALLOC(wt,g+1,double);

  copyarr(kmean, jmean, g) ; 
  copyarr(kmean, jmean+g, g) ; 
  copyarr(jwt, wt, g) ; 
  wt[g] = wt[0] ; 

  /**Calculate the mean of the jackknifed means*/
  mpr = asum(jmean, g) ; 
  mp = mpr / yg ; 

  gmean = mean ; 
 
  if (debug) printf("mp: %12.6f mean: %12.6f\n", mp, mean) ;

  vsp(jst, jmean, -mp, g) ;

  jst[g] = jst[0];

  rho = corr(jst, jst+1, g) ; 

  if (debug) printf("rho: %12.6f g:: %d\n", rho, g) ;  
  if (verbose) { 
   for (k=0; k<g; ++k) {
    printf("jst %d %12.6f %12.6f %12.6f %12.6f\n", k, jst[k], jmean[k], wt[k], mp);
   }
  }

  *prho = rho;

  if(rho < 0)
    {
     printf("(fourierjack rho negative!\n") ;   
     *prho = 0;
     weightjack(est,sig,mean,kmean,jwt,g) ; 
     free(jst);
     free(d);
     free(c);
     free(jp);
     free(wt);
     free(jmean) ; 
      return ; 
    }

  for(k=0;k<g;k++)
    {
      y = 1 + 2*rho*cos((2.0*pi*k)/yg);
      d[k] = 1.0/sqrt(y) ;
    }

  vzero(c, 2*g) ; 
  S = 0;
  for(i=0;i<g;i++)
    {
      for(k=0;k<g;k++)
	{
          y = (double) (k*i) ; 
	  c[i] += d[k]*cos((2.0*pi*y)/yg); 
	}
      c[i] = c[i]/yg;
      S += c[i];
    }
    y1 = asum(c, g) / yg ; 
    y2 = asum2(c, g) / yg ; 
//  printf("c moments %12.6f %12.6f\n", y1, y2 ) ;
//  printmat(c, 1, 40) ;  
    vsp(c, c, -y1, g) ; 
    y2 = asum2(c, g)  ; 
    vst(c, c, 1.0/sqrt(y2), g)  ;  
// c sums to zero and sum of squares is 1  

    y1 = asum(c, g)  ; 
    y2 = asum2(c, g) / yg ; 
//  printf("c adj S1 S2 %12.6f %12.6f\n", y1, y2 ) ;
//  printmat(c, 1, 40) ;  


  /*Using the periodicity of the cosine function*/
  copyarr(c, c+g, g) ; 
  copyarr(jst, jmean, g+1) ; 

  cS = 1.0/(sqrt(1+2.0*rho));

  jpmean = 0.0;
  for(i=0;i<g;i++)
    {
      jp[i] = vdot(c+i, jmean, g)  ;
    }
 
  // before tranform 

  jpmean = asum(jp, g)/yg;
/**
  y1 =  asum(jmean, g)/yg;
  y2 = asum2(jmean, g)/yg ; 
  y2 = sqrt(y2) ;  
  printf("A1 A2 %12.6f %12.6f\n", y1, y2) ;
  
  y1 = jpmean ;
  y2 = asum2(jp, g)/yg ; 
  y2 = sqrt(y2) ; 
  printf("B1 B2 %12.6f %12.6f\n", y1, y2) ;

  printmat(kmean, 1, 40) ;  printnl() ;
  printmat(jmean, 1, 40) ;  printnl() ;
*/

// overwrite jmean 


  copyarr(jp, jmean, g) ; 
  jmean[g] = jmean[0] ;

//  printmat(jmean, 1, 40) ;  printnl() ;

  ZALLOC(wk,2*g, double);
  ZALLOC(qq,2*g, double);
  y1 = asum(jmean, g) ; 
  y2 = (double) g ;  
  ymx = y1/y2  ; 
  copyarr(jmean, wk, g) ; 
  copyarr(wk, wk+g, g) ;  
  for (k=0; k<g; ++k) { 
   qq[k] = rho * sqrt(wt[k]*wt[k+1]) ; 
  } 
  y1 = asum(wk, g) ; 
  y2 = yg ; 
  vsp(wk, wk, -y1/y2, g) ; 
  wk[g] = wk[0] ; 
  ycx = corr(wk, wk+1, g) ; 
//  printf("zzchk %12.6f %12.6f %12.6f\n", ymx, y1/y2, ycx) ; 

  for (i=0; i<g; i++) { 
   y = 0; 
   for (k=0; k<g; ++k) { 
    l = k-i ; 
    if (l<0) l += g ;  
    y += c[k]*c[k]*wt[l] ; 
    y += c[k]*c[k+1]*qq[l] ; 
   }
   wk[i] = y ; 
  }
  printf("zz1\n") ; 
  printmat(wt, 1, 20) ;  
  printnl() ; 
  y = asum(wk, g) / yg ;  // mean 
  vst(wk, wk, 1.0/y, g) ; 
  copyarr(wk, wt, g) ; 
  free(wk) ; 
  free(qq) ; 

  vsp(jmean, jmean, mp, g) ; 
  weightjack(est, sig, gmean, jmean, wt, g) ;

//printf("S cS jpmean mean est %12.6f %12.6f %12.6f %12.5f %12.6f  rho:: %12.6f\n",S,cS,jpmean, gmean, *est, rho);
  
  free(jst);
  free(d);
  free(c);
  free(jp);
  free(wt);
  free(jmean) ; 
	     
}
