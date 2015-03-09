#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>

#define MAXSTR 512
#define MAXFF 20

char *iname = NULL ;
char *oname = NULL ;
double mean;

void readcommands(int argc, char **argv) ;
void callwtjack(char *iname, char *oname);

int main(int argc, char **argv)
{
  readcommands(argc, argv) ;
  callwtjack(iname, oname);
}

void readcommands(int argc, char **argv) 
{
  int i;

  while ((i = getopt (argc, argv, "i:o:m:")) != -1)
  {
    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'o':
	oname = strdup(optarg) ;
	break;

      case 'm':
         mean = atof(optarg) ;
	 break;

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }
}


void
callwtjack(char *iname, char *oname) 
{
  FILE *ifile, *ofile ;
  char line[MAXSTR];
  char *sx;
  char *spt[MAXFF];
  int nsplit, len, i, k;
  char c;
 
  openit(iname, &ifile, "r") ;
  openit(oname, &ofile, "w") ;

  double *jwt, *jmean;
  
  /* output variable */
  double est, sig;   // NB
  est =0; sig=0;

  /*input variables */
  len = numlines(iname);
  ZALLOC(jwt, len, double);     
  ZALLOC(jmean, len, double);
  k = 0;

 /* read input file and store data */
  while (fgets(line,MAXSTR,ifile) != NULL)  
  {
	nsplit = splitup(line, spt, MAXFF);
	sx = spt[0];
	c = sx[0];
	if (c == '#') {
		freeup(spt, nsplit);
		continue;
	}
	
	jwt[k] = atof(spt[1]);
	jmean[k] = atof(spt[2]);
	//printf("mean: %9.3f len: %9.3f\n", jmean[k], jwt[k]) ;
	k++;
	freeup(spt, nsplit);
  }  
  len = k ;  // better style  who knows how numlines handles commas
 fclose(ifile) ;
 // printf("mean: %9.3f len: %d\n", mean, len) ;


 /*call weightjack */
	weightjack(&est, &sig, mean, jmean, jwt, len);
	fprintf(ofile,"%9.3f", est);	 	// d format ??
	fprintf(ofile,"%9.3f", sig);	 	
	fprintf(ofile,"\n");
	free(jmean);
	free(jwt);
 fclose(ofile) ;

 }
