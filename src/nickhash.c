#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>
#include "mcmcpars.h"
#include "globals.h"
#include "admutils.h"

#define WVERSION   "100" 

extern int verbose  ;
int qtmode = NO ;
int dolower = NO ;

char *iname = NULL ;

void readcommands(int argc, char **argv) ;  


int main(int argc, char **argv)
{

  char **names ; 
  int n, z ; 
  char ss[10] ;

// printf("hello world!\n") ;

  readcommands(argc, argv) ;
  if (iname == NULL) {
    printf("Usage: nickhash -i <yourfile> [-l]\n") ;
    return -1 ; 
  }

  n = numlines(iname) ; 
  ZALLOC(names, n, char *) ;
  n = getlist(iname, names) ; 

  z = hasharr(names, n) ; 
  sprintf(ss, "%08X", z) ;
  if (dolower) mklower(ss) ;
  printf("%s", ss) ;
  
  printnl() ;
 
  return 0 ;

}

void readcommands(int argc, char **argv) 

{
  int i ;

  while ((i = getopt (argc, argv, "i:lvV")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
        break ; 

      case 'V':
        verbose = YES ;
        break ; 

      case 'l':
        dolower = YES ;
        break ; 


     }

  }
 
  return ;

}
