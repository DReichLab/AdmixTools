#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <nicklib.h>
#include <getpars.h>

#define MAXSTR 512
#define MAXFF 20

char *dir = "/home/np29/broaddata/tables/hg17";
char *iname = NULL;
char *oname = NULL;

int verbose = NO;
char *grabp (int argc, char **argv);

int
main (int argc, char **argv)
{

  char *xx;

  xx = grabp (argc, argv);
  if (xx == NULL)
    xx = strdup ("NOTFOUND");
  printf ("%s\n", xx);
  if (strlen (xx) > 0)
    return 0;

  return -1;


}

char *
grabp (int argc, char **argv)
{
  int i;
  phandle *ph;
  char str[512];
  int n, kode;
  char *parname = NULL, *xname = NULL;
  static char *xval = NULL;

  while ((i = getopt (argc, argv, "p:x:V")) != -1) {

    switch (i) {

    case 'p':
      parname = strdup (optarg);
      break;

    case 'x':
      xname = strdup (optarg);
      break;

    case 'V':
      verbose = YES;
      break;

    case '?':
      printf ("Usage: bad params.... \n");
      fatalx ("bad params\n");
    }
  }

  if (parname == NULL)
    return NULL;
  if (xname == NULL)
    return NULL;

  ph = openpars (parname);
  dostrsub (ph);

  getstring (ph, xname, &xval);
  return xval;

}
