#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <libgen.h>

#include <nicklib.h>
#include <getpars.h>

#define MAXSTR 512
#define MAXFF 20

char *dir = "/home/np29/broaddata/tables/hg17";
char *iname = NULL;
char *oname = NULL;

int verbose = NO;
char *grabp (int argc, char **argv);

int usage (char *prog, int exval);

int usage (char *prog, int exval)
{

  (void)fprintf(stderr, "Usage: %s [options] -x <nam> -p <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters <file> .\n");
  (void)fprintf(stderr, "   -x <nam>    ... extract value of parameter <nam> from param file.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");

  exit(exval);
}

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

  if (argc == 1) { usage(basename(argv[0]), 1); }

  while ((i = getopt (argc, argv, "hp:x:V")) != -1) {

    switch (i) {

    case 'h':
      usage(basename(argv[0]), 0);

    case 'p':
      parname = strdup (optarg);
      break;

    case 'x':
      xname = strdup (optarg);
      break;

    case 'V':
      verbose = YES;
      break;

    default:
      usage(basename(argv[0]), 1);
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
