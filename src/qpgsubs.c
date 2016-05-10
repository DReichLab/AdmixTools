#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>
#include "regsubs.h"
#include "qpsubs.h"
#include "admutils.h"

#define WVERSION "2000"

extern int verbose;
extern int isinit;
extern int hires;
static int totpop = 0;
static int *ispath = NULL;
static int ispathlen = 0;
// ispath [x*numvertex=y] = 1 iff path from x -> y

#define HUGEDIS (1 << 29)


typedef struct
{
  struct NODE *vlist;
  int root;
} GRAPH;

typedef struct
{
  int gnode;
  int glabel;
  char *name;
  char *label;
  int isroot;
  int isleaf;
  int distance;
  int isadmix;
  int isfixed;                  //  1 => no reestimate
  int numadmix;                 //index in adlist
  struct NODE *left;
  struct NODE *right;
  struct NODE *parent;
  struct EDGE *eleft;
  struct EDGE *eright;
  struct EDGE *eparent;
  double time;
  int popsize;
  int numwind;
  int windex[MAXW];
  double wmix[MAXW];
  int supernum;
  int isdead;
  int eglistnum;
} NODE;

typedef struct
{
  char *name;
  struct NODE *up;
  struct NODE *down;
  double val;
  int edgenum;
  int isfixed;
  int iszero;
  int isleft;
} EDGE;

static NODE *vlist = NULL;
static EDGE *elist = NULL;
static char **eglist = NULL;

static int *egnum;              // vlist number correspnding to pop
static NODE **adlist = NULL;
static NODE **xadlist = NULL;   // variable admix
static char **forcenames;       // edgenames forced to zero
static int nforce = 0;

static int numedge = -1, numvertex = -1, numadmix, numpops;
int ncall;

void initelist (EDGE * elist, int n);
void initvlist (NODE * vlist, int n);
int readit (char *cname);
int qreadit (char *cname);
int vindex (char *vname, NODE * vlist, int n);
void addwts (double *vv, NODE * node, double weight);
void addadmix (char *vertname, char **avnames, double *ww, int nt);
int getadwts (char **spt, char **ssx, double *ww, int nsplit, int *n);
int getextlabel (NODE * node, char *sss);
void ckline (char *line, int ns, int n);

EDGE *xedge (NODE * n1, NODE * n2);
void fixnode (NODE * node);
void fixedge (EDGE * edge);
NODE *root ();
NODE *vert (char *nodename);
void setdistances (NODE * node);
void cleardis ();
void setdis (NODE * node, int dis);
int setxx (NODE * node, NODE ** xx);
void reroot (char *nodename);
void settime ();
void setvnum (int *vnum, int *list, int n);
void pedge (FILE * fff, NODE * anode, NODE * bnode, double val, int mode);
int vertexnum (char *vertname);
int edgenum (char *edgename);

void
getvnames (char **vnames)
{
  int k;
  NODE *node;

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    vnames[k] = strdup (node->name);
  }
}

void
findvname (char *ss, int n)
{
  char **vnames;
  char sss[20];
  int x;

  ZALLOC (vnames, numvertex, char *);
  getvnames (vnames);

  x = 1;

  for (;;) {
    sprintf (sss, "V%dX%d", n, x);
    if (indxstring (vnames, numvertex, sss) < 0)
      break;
    ++x;
  }

  strcpy (ss, sss);

  freeup (vnames, numvertex);
  free (vnames);

  return;
}

void
findename (char *ss, int n)
{
  char **enames;
  char sss[20];
  int x;

  ZALLOC (enames, numedge, char *);
  getenames (enames);

  x = 1;

  for (;;) {
    sprintf (sss, "e%dx%d", n, x);
    if (indxstring (enames, numedge, sss) < 0)
      break;
    ++x;
  }

  strcpy (ss, sss);


  freeup (enames, numedge);
  free (enames);

  return;
}

int
getedgelock (int *lock, double *vals)
{
  int k, x = 0;
  EDGE *edge;

  vzero (vals, numedge);
  ivclear (lock, -1, numedge);
  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    if (edge->isfixed == NO)
      continue;
    vals[x] = edge->val;
    lock[x] = k;
    ++x;
  }
  return x;
}


void
getenames (char **enames)
{
  int k;
  EDGE *edge;

  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    enames[k] = strdup (edge->name);
  }
}

void
getezero (int *zpat)
{
  int k;
  EDGE *edge;

  ivzero (zpat, numedge);
  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    if (edge->iszero)
      zpat[k] = 1;
  }

}

int
grdof ()
{

  int t, k, l;
  NODE *node;

  if (numedge <= 0)
    fatalx ("no graph\n");
  t = numedge;
  for (k = 0; k < numadmix; ++k) {
    node = adlist[k];
    l = node->numwind;
    if (l == 0)
      continue;
    t += (l - 1);
  }

  return t;

}

void
getmixstr (int k, char *sss)
{
  NODE *node, *nodex;
  int j, l, t;

  node = xadlist[k];

  sss[0] = CNULL;
  getextlabel (node, sss);
  strcat (sss, " ");
  strcat (sss, node->name);
  l = node->numwind;
  for (j = 0; j < l; j++) {
    t = node->windex[j];
    nodex = &vlist[t];
    strcat (sss, " ");
    strcat (sss, nodex->name);
  }
}

int
getrootlabel (char *sss)
{
  NODE *node;

  sss[0] = CNULL;
  node = root ();
  return getextlabel (node, sss);

}

int
getextlabel (NODE * node, char *sss)
{

  if (node == NULL)
    return;
  if (node->label != NULL) {
    strcat (sss, node->label);
    return 1;
  }
  if ((node->left) != NULL) {
    strcat (sss, ":");
    getextlabel ((NODE *) node->left, sss);
    return 0;
  }
  if ((node->right) != NULL) {
    strcat (sss, ":");
    getextlabel ((NODE *) node->right, sss);
    return 0;
  }
}

void
getvmind (int ind, int *adind, int *wind)
{
  int t, k, j, l, z, lo, hi;
  NODE *node;

  t = z = 0;

  for (k = 0; k < numadmix; ++k) {
    node = adlist[k];
    l = node->numwind;
    if (l == 0)
      continue;
    lo = z;
    hi = z + l - 1;
    if (hi >= ind) {
      *adind = t;
      *wind = ind - lo;
      return;
    }
    z += l;
    ++t;
  }
  fatalx ("bad getvmind\n");
}

int
getnumadmix ()
{
  if (numadmix < 0)
    fatalx ("graph not initialized\n");
  return numadmix;
}

int
getnumvertex ()
{
  if (numvertex < 0)
    fatalx ("graph not initialized\n");
  return numvertex;
}

int
getnumedge ()
{
  if (numedge < 0)
    fatalx ("graph not initialized\n");
  return numedge;
}

void
getewts (double *ewts)
{
  int k;
  EDGE *edge;
  NODE *node;
  char *s1, *s2;
  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    ewts[k] = edge->val;
    if (isnan (ewts[k])) {
      node = (NODE *) edge->up;
      s1 = node->name;
      node = (NODE *) edge->down;
      s2 = node->name;
      printf ("bad edge: %s %s %s\n", edge->name, s1, s2);
    }
  }
}

void
putewts (double *ewts)
{
  int k;
  EDGE *edge;
  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    edge->val = MAX (ewts[k], 0.0);
  }
}

void
getvwts (double *pwts, int *nrows, int *nedge)
// pwts preallocated like getpwts but for all vertices
{
  double *vv;
  int i, k;
  NODE *node;
  double *wts;

  wts = pwts;
  ZALLOC (vv, numedge, double);
  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    ncall = 0;
    vzero (vv, numedge);
    addwts (vv, node, 1.0);
    copyarr (vv, wts, numedge);
    wts += numedge;
  }
  free (vv);
  *nrows = numvertex;
  *nedge = numedge;
}

int
getpopinfo (char *name, int *popsize, int num)
{
  NODE *node;
  int k;

  k = egnum[num];
  node = &vlist[k];
  *popsize = node->popsize;
  strcpy (name, node->name);

  return k;

}

void
getpwtsx (double *pwts, int *nrows, int *nedge)
// pwts preallocated
// like getpwts but no special treatment of base
{
  double *vv;
  int i, k;
  NODE *node;
  double *wts;

  wts = pwts;
  ZALLOC (vv, numedge, double);
  for (i = 0; i < numpops; ++i) {
    k = egnum[i];
    node = &vlist[k];
    ncall = 0;
    vzero (vv, numedge);
    addwts (vv, node, 1.0);
    copyarr (vv, wts, numedge);
    wts += numedge;
  }
  free (vv);
  *nrows = numpops;
  *nedge = numedge;
}

void
getpwts (double *pwts, int *nrows, int *nedge)
// pwts preallocated
{
  double *vv;
  int i, k;
  NODE *node;
  double *wts;

  wts = pwts;
  ZALLOC (vv, numedge, double);
  for (i = 1; i < numpops; ++i) {
// base population is zero
    k = egnum[i];
    node = &vlist[k];
    ncall = 0;
    vzero (vv, numedge);
    addwts (vv, node, 1.0);
    copyarr (vv, wts, numedge);
    wts += numedge;
  }
  free (vv);
  *nrows = numpops - 1;
  *nedge = numedge;
}

void
getpnwts (double *pwts, int *nrows, int *nedge)
// pwts preallocated.  Like pnwts but for all pops
{
  double *vv;
  int i, k;
  NODE *node;
  double *wts;

  wts = pwts;
  ZALLOC (vv, numedge, double);
  for (i = 0; i < numvertex; ++i) {
// base population is zero
    node = &vlist[i];
    ncall = 0;
    vzero (vv, numedge);
    addwts (vv, node, 1.0);
    copyarr (vv, wts, numedge);
    wts += numedge;
  }
  free (vv);
  *nrows = numvertex - 1;
  *nedge = numedge;
}

void
addwts (double *vv, NODE * node, double weight)
// recursively fill in generalized path weights
{
  NODE *tnode;
  EDGE *edge;
  int k, t, vind;
  double wt;

  ++ncall;
  if (ncall > 100) {
    fatalx ("looping\n");
  }

  if (node->isroot)
    return;
  t = node->numwind;
  if (t > 0) {
    for (k = 0; k < t; ++k) {
      wt = node->wmix[k];
      vind = node->windex[k];
      tnode = &vlist[vind];
      addwts (vv, tnode, weight * wt);
    }
    return;
  }
  edge = (EDGE *) node->eparent;
  if (edge == NULL) {
    fatalx ("node: %s has no parent\n", node->name);
  }
  t = edge->edgenum;
  vv[t] += weight;
  tnode = (NODE *) edge->up;
  addwts (vv, tnode, weight);
}

void
getgmix (double **vmix, int *lmix, int *nmix)
// vmix is [100][MAXW] and preallocated
// extract mixing coefficients  
{
  int t, k, l;
  NODE *node;

  t = 0;
  clear2D (&vmix, numadmix, MAXW, 0.0);
  for (k = 0; k < numadmix; ++k) {
    node = adlist[k];
    l = node->numwind;
    if (l == 0)
      continue;
    if (node->isfixed)
      continue;
    xadlist[t] = node;
    copyarr (node->wmix, vmix[t], l);
    lmix[t] = l;
    ++t;
  }
  *nmix = t;
}

void
putgmix (double **vmix)
// shape of mixers already in  nodes of vlist
// weights are forced to sum to 1
{
  int t, k, l;
  NODE *node;

  t = 0;
  for (k = 0; k < numadmix; ++k) {
    node = adlist[k];
    l = node->numwind;
    if (l == 0)
      continue;
    if (node->isfixed)
      continue;
    vzero (node->wmix, MAXW);
    vabs (node->wmix, vmix[t], l);
    bal1 (node->wmix, l);
    ++t;
  }
}

void
destroyg ()
{

  if (vlist != NULL)
    free (vlist);
  if (elist != NULL)
    free (elist);
  if (adlist != NULL)
    free (adlist);
  if (xadlist != NULL)
    free (xadlist);

  if (numpops != 0) {
    freeup (eglist, numpops);
    free (eglist);
  }

  if (ispath != NULL) {
    free (ispath);
    ispath = NULL;
  }

  ispathlen = 0;
  numedge = numvertex = numadmix = numpops = 0;

}

void
setispath ()
{
  int *w1, *w2;
  int n = numvertex;
  int nn, a, b;

  nn = n * n;
  ZALLOC (w1, nn, int);
  ZALLOC (w2, nn, int);

  setincidence (w1);
  ivclip (w1, w1, 0, 1, nn);
  a = intsum (w1, nn);
  if (a == 0)
    fatalx ("bug\n");
  for (;;) {
    imulmat (w2, w1, w1, n, n, n);
    ivvp (w2, w2, w1, nn);
    ivclip (w1, w2, 0, 1, nn);
    b = intsum (w1, nn);
    if (a == b)
      break;
    if (a < b)
      fatalx ("badbug A\n");
    if (b > nn)
      fatalx ("badbug B\n");
    a = b;
  }
  if (ispath != NULL)
    free (ispath);
  ZALLOC (ispath, nn, int);
  copyiarr (w1, ispath, nn);
  free (w1);
  free (w2);
}

void
setincidence (int *x)
// not clipped
{
  int n = numvertex;
  int a, u, v;
  int k, t;
  NODE *node;
  EDGE *edge;

  ivzero (x, n * n);
  for (a = 0; a < numedge; a++) {
    edge = &elist[a];
    node = (NODE *) edge->up;
    u = node->gnode;
    node = (NODE *) edge->down;
    v = node->gnode;
    x[u * n + v] += 1;
  }
  for (a = 0; a < numvertex; a++) {
    node = &vlist[a];
    t = node->numwind;
    v = a;
    for (k = 0; k < t; ++k) {
      u = node->windex[k];
      x[u * n + v] += 1;
    }
    u = v = a;
    x[u * n + v] = 1;
  }
}

void
trivgraph (char *rootname)
{
  int maxnodes = MAXG;
  int numeg, k;
  NODE *node;
  EDGE *ep;

  destroyg ();

  ZALLOC (vlist, maxnodes, NODE);
  ZALLOC (elist, maxnodes, EDGE);
  ZALLOC (adlist, maxnodes, NODE *);
  ZALLOC (xadlist, maxnodes, NODE *);
  ZALLOC (eglist, maxnodes, char *);

  initvlist (vlist, maxnodes);
  initelist (elist, maxnodes);

  addvertex (rootname);
}

void
addadmix (char *vertname, char **avnames, double *ww, int nt)
{

  NODE *node;
  int t, k, tbase;

  for (k = 0; k < nt; ++k) {
    addvertex (avnames[k]);
  }
  addvertex (vertname);
  tbase = t = vertexnum (vertname);
  node = &vlist[t];
  for (k = 0; k < nt; ++k) {
    t = vertexnum (avnames[k]);
    if (t == tbase)
      fatalx ("bad admix (recursion) %s\n", vertname);
    node->windex[k] = t;
    ++(node->numwind);
  }
  copyarr (ww, node->wmix, nt);
  node->isadmix = YES;
  adlist[numadmix] = node;
  node->numadmix = numadmix;
  ++numadmix;
}

int
gsimplify (int n)
{
  int nkill = 0, i, k, l, t;
  NODE *node, *nparent, *tnode;
  int *ecount;
  char tmpqq[20];
  double val, vparent, vchild;
  EDGE *cedge, *pedge;

  if (n > 100)
    fatalx ("looping\n");
  for (i = 0; i < numvertex; ++i) {
    node = &vlist[i];
    if (node->isdead)
      continue;
    l = node->numwind;
    if (l == 0)
      continue;
    for (k = 0; k < l; ++k) {
      t = node->windex[k];
      tnode = &vlist[t];
      tnode->isleaf = NO;
    }
  }

  for (i = 0; i < numvertex; ++i) {
    node = &vlist[i];
    if (!node->isleaf)
      continue;
    if (node->label != NULL)
      continue;
    if (node->isdead)
      continue;
    ++nkill;
// printf("##deleting %s\n", node -> name) ;
    node->isdead = YES;
    nparent = (NODE *) node->parent;
    if (nparent == NULL)
      continue;
    node->isleaf = NO;
    node->numwind = 0;

    if ((NODE *) nparent->left == node) {
      nparent->left = NULL;
      nparent->eleft = NULL;
      if (nparent->right == NULL)
        nparent->isleaf = YES;
    }

    if ((NODE *) nparent->right == node) {
      nparent->right = NULL;
      nparent->eright = NULL;
      if (nparent->left == NULL)
        nparent->isleaf = YES;
    }
    node->left = node->right = NULL;
    node->eleft = node->eright = NULL;
  }
  if (nkill > 0)
    return gsimplify (n + nkill);

  ZALLOC (ecount, numvertex, int);
  for (i = 0; i < numvertex; ++i) {
    node = &vlist[i];
    if (node->isdead)
      continue;
    l = node->numwind;
    for (k = 0; k < l; ++k) {
      t = node->windex[k];
      ecount[t] += 10;
    }
    if (node->right != NULL)
      ++ecount[i];
    if (node->left != NULL)
      ++ecount[i];
  }
// kill edges with one outedge
  for (i = 0; i < numvertex; ++i) {
    node = &vlist[i];
    if (node->isroot)
      continue;
    if (node->isdead)
      continue;
    if (node->numwind > 0)
      continue;
    if (ecount[i] != 1)
      continue;
    ++nkill;
//  printf("##deleting %s\n", node -> name) ;
    node->isdead = YES;
    pedge = (EDGE *) node->eparent;
    if (node->left != NULL) {
      tnode = (NODE *) node->left;
      cedge = (EDGE *) node->eleft;
    }
    if (node->right != NULL) {
      tnode = (NODE *) node->right;
      cedge = (EDGE *) node->eright;
    }

    pedge->val += cedge->val;

    nparent = (NODE *) node->parent;
    tnode->parent = (struct NODE *) nparent;
    if ((NODE *) nparent->left == node) {
      nparent->left = (struct NODE *) tnode;
    }
    if ((NODE *) nparent->right == node) {
      nparent->right = (struct NODE *) tnode;
    }
    node->left = node->right = NULL;
    node->eleft = node->eright = NULL;
  }
  free (ecount);
  sprintf (tmpqq, "graphtt:%d", getpid);
  dumpgraph (tmpqq);
  loadgraph (tmpqq, NULL);
  unlink (tmpqq);
  return n + nkill;
}

int
dellabel (char *label)
{
  NODE *node;
  int t, i;
  char **tlist;

  t = findlabel (label);
  if (t < 0)
    return -99;
  node = &vlist[t];
  if (node->isroot)
    return -1;

  freestring (&node->label);
  t = indxindex (eglist, numpops, label);
  if (t < 0)
    fatalx ("bugdel\n");

  ZALLOC (tlist, numpops - 1, char *);
  copystringsd (eglist, tlist, numpops, t);
  freeup (eglist, numpops);
  copystrings (tlist, eglist, numpops - 1);
  freeup (tlist, numpops - 1);
  --numpops;

  return numpops;

}

void
addlabel (char *label, char *vertname)
{
  NODE *node;
  int t, i;

  addvertex (vertname);
  t = vertexnum (vertname);
  node = &vlist[t];
  node->label = strdup (label);
  eglist[numpops] = strdup (label);
//  printf("zz %d %s\n", numpops, label) ;
  ++numpops;
}

void
addvertex (char *vertname)
{
  NODE *node;
  int t;

  if (vertname == NULL)
    fatalx ("bad addvertex call\n");
  t = vertexnum (vertname);
  if (t >= 0)
    return;
  node = &vlist[numvertex];
  node->name = strdup (vertname);
  ++numvertex;
}

void
addedge (char *ename, char *x1, char *x2)
{
  EDGE *edge;
  NODE *node1, *node2, *node, *nl, *nr;
  int t, s1, s2;

  addvertex (x1);
  addvertex (x2);

  t = edgenum (ename);
  if (t >= 0)
    fatalx ("edge %s exists!\n", ename);
  edge = &elist[numedge];
  edge->name = strdup (ename);
  // printf("addedge: %s %d\n", edge -> name, numedge) ;  fflush(stdout) ;
  ++numedge;
  s1 = vertexnum (x1);
  s2 = vertexnum (x2);
  node2 = &vlist[s2];
  node = node1 = &vlist[s1];

  nl = (NODE *) node1->left;
  nr = (NODE *) node1->right;

  if ((node->eleft) == NULL) {
    node->eleft = (struct EDGE *) edge;
    edge->isleft = YES;
    node = &vlist[s2];
    node->eparent = (struct EDGE *) edge;
    edge->up = (struct NODE *) node1;
    edge->down = (struct NODE *) node2;
    node1->left = (struct NODE *) node2;
    node2->parent = (struct NODE *) node1;
    node1->isleaf = NO;
    return;
  }
  if ((node->eright) != NULL)
    fatalx ("> 2 edges at %s\n", node->name);
  node->eright = (struct EDGE *) edge;
  edge->isleft = NO;
  node = &vlist[s2];
  node->eparent = (struct EDGE *) edge;
  edge->up = (struct NODE *) node1;
  edge->down = (struct NODE *) node2;
  node1->right = (struct NODE *) node2;
  node2->parent = (struct NODE *) node1;
  node1->isleaf = NO;
}

void
freegraph ()
{
  int i;
  NODE *node;
  EDGE *edge;

  if (numvertex <= 0)
    return;
  for (i = 0; i < numvertex; ++i) {
    node = &vlist[i];
    freestring (&node->name);
    freestring (&node->label);
  }
  free (vlist);
  for (i = 0; i < numedge; ++i) {
    edge = &elist[i];
    freestring (&edge->name);
    freestring (&forcenames[i]);
  }
  for (i = 0; i < numpops; ++i) {
    freestring (&eglist[i]);
  }
  if (eglist != NULL)
    free (eglist);
  eglist = NULL;

  free (elist);
  free (adlist);
  free (eglist);
  free (forcenames);
  free (egnum);

  vlist = NULL;
  elist = NULL;

  numedge = numvertex = numadmix = numpops = -1;
}

int
loadeglist (char ***pxeglist, int xnumeg)
{

  static char **xeglist;
  char **teglist;
  int numeg = numpops, k;
  char *tlist;

  teglist = *pxeglist;


  if (teglist != NULL) {
    freeup (teglist, xnumeg);
    free (teglist);
  }

  ZALLOC (xeglist, numpops, char *);

/**
  for (k=0; k<numpops; ++k) {
   printf("zz2 %s\n", eglist[k]) ;
  }
  fflush(stdout) ;
*/

  for (k = 0; k < numpops; ++k) {
    xeglist[k] = strdup (eglist[k]);
  }

  *pxeglist = xeglist;
  return numpops;

}

int
loadgraph (char *rname, char ***peglist)
{
  int maxnodes = MAXG;
  int numeg, k;
  NODE *node;
  EDGE *ep;

  freegraph ();

  numedge = numvertex = numadmix = numpops = 0;

  ZALLOC (egnum, maxnodes, int);
  ZALLOC (vlist, maxnodes, NODE);
  ZALLOC (elist, maxnodes, EDGE);
  ZALLOC (adlist, maxnodes, NODE *);
  ZALLOC (xadlist, maxnodes, NODE *);
  ZALLOC (forcenames, maxnodes, char *);
  initvlist (vlist, maxnodes);
  initelist (elist, maxnodes);
  numvertex = readit (rname);

  node = root ();
  node->isroot = YES;
  // printf("numvertex: %d  root: %s\n", numvertex, node -> name) ;
  numeg = 0;
  ZALLOC (eglist, MAXG, char *);
  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    if (node->label == NULL)
      continue;
    eglist[numeg] = strdup (node->label);
    egnum[numeg] = k;
    node->eglistnum = numeg;
    ++numeg;
  }

  numpops = numeg;

  if (peglist != NULL) {
    ZALLOC (*peglist, numeg, char *);
    copystrings (eglist, *peglist, numeg);
  }

/**
  for (k=0; k<numvertex; ++k)  {  
   node = &vlist[k] ;
   printf("zzp: %s ", node -> name) ;
   ep = (EDGE *) node -> eparent ;  
   if (ep != NULL) printf(" %s", ep ->name) ;
   printnl() ;
  }
*/

  return numeg;
}

void
initelist (EDGE * elist, int n)
{
  EDGE *edge;
  int k;

  for (k = 0; k < n; ++k) {
    edge = &elist[k];
    edge->up = edge->down = NULL;
    edge->edgenum = k;
    edge->name = NULL;
    edge->val = 0.0;
    edge->iszero = NO;
    edge->isleft = NO;
    edge->isfixed = NO;
  }
}

void
initvlist (NODE * vlist, int n)
{
  NODE *node;
  int k;

  for (k = 0; k < n; ++k) {
    node = &vlist[k];
    node->gnode = k;
    node->name = node->label = NULL;
    node->isroot = NO;
    node->isleaf = YES;
    node->isadmix = NO;
    node->isfixed = NO;
    node->numadmix = -1;
    node->numwind = 0;
    node->popsize = 0;
    node->time = 0.0;
    node->supernum = -1;
    node->eglistnum = -1;
    node->isdead = NO;

    ivclear (node->windex, -1, MAXW);
    node->left = node->right = node->parent = NULL;
  }
}

int
readit (char *cname)
// main input routine for graph
{
#define MAXSTR 128
#define MAXFF  10
  FILE *fff;
  char line[MAXSTR + 1], c;
  char *spt[MAXFF], *sx, *s1, *s2;
  int nsplit, n = 0;
  int kret, t, v1, v2, w, k;
  int okline;
  NODE *node, *tnode, *node2, *node1;
  EDGE *edge;
  double ww[MAXW];
  char *ssx[MAXW];
  int nt;
  double val;
  int num = 0;


  openit (cname, &fff, "r");
  line[MAXSTR] = '\0';
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    ++num;

    kret = strcmp (sx, "root");
    if ((kret == 0) && (num == 1)) {
      fclose (fff);
      return qreadit (cname);
    }

    okline = NO;
    kret = strcmp (sx, "vertex");
    if (kret == 0) {
      node = &vlist[n];
      node->name = strdup (spt[1]);
      ++n;
      okline = YES;
    }
    kret = strcmp (sx, "label");
    if (kret == 0) {
      sx = spt[1];
      t = vindex (sx, vlist, n);
      if (t < 0)
        fatalx ("bad label: %s\n", sx);
      sx = spt[2];
      tnode = &vlist[t];
      tnode->label = strdup (sx);
      if (nsplit >= 4) {
        sx = spt[3];
        tnode->popsize = atoi (sx);
      }
      okline = YES;
    }
    kret = strcmp (sx, "lock");
    if (kret == 0) {
      sx = spt[1];
      t = vindex (sx, vlist, n);
      if (t < 0)
        fatalx ("bad lock: %s\n", sx);
      tnode = &vlist[t];
      tnode->isfixed = YES;
      okline = YES;
    }
    kret = strcmp (sx, "ledge");
    if (kret == 0) {
      if (nsplit<4) fatalx("no target node: %s\n", line) ;
      edge = &elist[numedge];
      sx = spt[1];
      s1 = spt[2];
      s2 = spt[3];
      edge->name = strdup (sx);
      if (nsplit == 6) {
        sx = spt[5];
        t = strcmp (sx, "lock");
        if (t != 0)
          fatalx ("bad line %s\n", line);
        edge->isfixed = YES;
      }
      if (nsplit >= 5) {
        sx = spt[4];
        t = strcmp (sx, "zero");
        if (t == 0)
          edge->iszero = YES;
        if (!isalpha (sx[0]))
          edge->val = atof (sx);
      }
      v1 = vindex (s1, vlist, n);
      v2 = vindex (s2, vlist, n);
      if (v1 < 0)
        fatalx ("bad edge: %s\n", line);
      if (v2 < 0)
        fatalx ("bad edge: %s\n", line);
      node1 = &vlist[v1];
      if (node1->left != NULL)
        fatalx ("two left edges! %s\n", node1->name);
      node2 = &vlist[v2];
      node = &vlist[v1];
      node->eleft = (struct EDGE *) edge;
      edge->isleft = YES;
      node = &vlist[v2];
      node->eparent = (struct EDGE *) edge;
      edge->up = (struct NODE *) node1;
      edge->down = (struct NODE *) node2;
      node1->left = (struct NODE *) node2;
      node2->parent = (struct NODE *) node1;
      node1->isleaf = NO;
      ++numedge;
      okline = YES;
      edge->val = MAX (edge->val, 0.0);
    }
    kret = strcmp (sx, "redge");
    if (kret == 0) {
      if (nsplit<4) fatalx("no target node: %s\n", line) ;
      edge = &elist[numedge];
      if (nsplit == 6) {
        sx = spt[5];
        t = strcmp (sx, "lock");
        if (t != 0)
          fatalx ("bad line %s\n", line);
        edge->isfixed = YES;
      }
      if (nsplit >= 5) {
        sx = spt[4];
        t = strcmp (sx, "zero");
        if (t == 0)
          edge->iszero = YES;
        if (!isalpha (sx[0]))
          edge->val = atof (sx);
      }
      sx = spt[1];
      s1 = spt[2];
      s2 = spt[3];
      edge->name = strdup (sx);
      v1 = vindex (s1, vlist, n);
      v2 = vindex (s2, vlist, n);
      if (v1 < 0)
        fatalx ("bad edge: %s\n", line);
      if (v2 < 0)
        fatalx ("bad edge: %s\n", line);
      node1 = &vlist[v1];
      if (node1->right != NULL)
        fatalx ("two right edges! %s\n", node1->name);
      node2 = &vlist[v2];
      node = &vlist[v1];
      node->eright = (struct EDGE *) edge;
      node = &vlist[v2];
      node->eparent = (struct EDGE *) edge;
      edge->up = (struct NODE *) node1;
      edge->down = (struct NODE *) node2;
      node1->right = (struct NODE *) node2;
      node2->parent = (struct NODE *) node1;
      node1->isleaf = NO;
      ++numedge;
      okline = YES;
      edge->val = MAX (edge->val, 0.0);
    }
    kret = strcmp (sx, "admix");
    if (kret == 0) {
      s1 = spt[1];
      v1 = vindex (s1, vlist, n);
      if (v1 < 0)
        fatalx ("bad admix: %s\n", line);
      node = &vlist[v1];
      sx = spt[nsplit - 1];
      t = strcmp (sx, "lock");
      if (t == 0) {
        node->isfixed = 1;
        getadwts (spt, ssx, ww, nsplit - 1, &nt);
      }
      else
        getadwts (spt, ssx, ww, nsplit, &nt);
      w = 0;
      for (k = 0; k < nt; ++k) {
        s2 = ssx[k];
        v2 = vindex (s2, vlist, n);
        if (v2 < 0)
          fatalx ("bad admix: %s\n", line);
        if (v1 == v2)
          fatalx ("bad admix (recursion) %s\n", line);
        node->windex[k] = v2;
        ++(node->numwind);
      }
      copyarr (ww, node->wmix, nt);
      node->isadmix = YES;
      adlist[numadmix] = node;
      ++numadmix;
      okline = YES;
    }
    if (okline == NO) {
      printf ("*** warning bad line:\n %s\n", line);
    }
    freeup (spt, nsplit);
    continue;
  }


  fclose (fff);
  return n;

}

int
getadwts (char **spt, char **ssx, double *ww, int nsplit, int *n)
{

  int nt, nw, k;
  char *sx;

  nt = 0;
  nw = 0;
  for (k = 2; k < nsplit; ++k) {
    sx = spt[k];
    if (isalpha (sx[0])) {
      ssx[nt] = sx;
      ++nt;
    }
    else {
      ww[nw] = atof (sx);
      ++nw;
    }
  }
  if (nw == 0) {
    setsimp (ww, nt);
    *n = nt;
    return;
  }
  if (nw != nt) {
    fatalx ("input wts for node %s not set correctly\n", spt[1]);
  }

  bal1 (ww, nt);                // balance weights
  vclip (ww, ww, .001, 1.0, nt);
  bal1 (ww, nt);                // balance weights
  isinit = YES;
  *n = nt;
  return;
}

int
popvindex (char *popname)
{
  int k, t;
  char *lab;

  for (k = 0; k < numvertex; ++k) {
    lab = vlist[k].label;
    if (lab == NULL)
      continue;

    t = strcmp (popname, lab);
    if (t == 0)
      return k;
  }
  return -1;
}

int
vind (char *vname)
{
  return vindex (vname, vlist, numvertex);
}


int
vindex (char *vname, NODE * vlist, int n)
{
  int k, t;

  for (k = 0; k < n; ++k) {
    t = strcmp (vname, vlist[k].name);
    if (t == 0)
      return k;
  }
  return -1;
}

void
setrand (double *ww, int n)
// fill array with uniform  
{
  int k;

  for (k = 0; k < n; ++k) {
    ww[k] = DRAND ();
  }
}

void
setsimp (double *ww, int n)
{
  double *pp;
  ZALLOC (pp, n, double);
  vclear (pp, 1.0, n);

  randirichlet (ww, pp, n);

  free (pp);

}

void
setpopsizes (int *sizes, char **eglist, int numeg)
{

  NODE *node;
  int k, t;

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    if (node->label == NULL)
      continue;
    t = indxindex (eglist, numeg, node->label);
    if (t < 0)
      continue;
    node->popsize = sizes[t];
    totpop += sizes[t];
  }
}

void
dumpdotgraph (char *graphdotname)
// in dot format 
{

  FILE *fff;
  NODE *node, *xnode, *xroot;
  EDGE *edge;
  int k, j, t, vind, kk;
  int *dd, *ind;
  char sform[10];
  double val, vmax;
  char *sss = NULL;

  if (graphdotname == NULL)
    return;
  printf ("graphdotname:  %s\n", graphdotname);
  fflush (stdout);
  openit (graphdotname, &fff, "w");

/** 
 we process vertices from the root 
*/
  xroot = root ();
  setdistances (xroot);
  ZALLOC (dd, numvertex, int);
  ZALLOC (ind, numvertex, int);
  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    dd[k] = node->distance;
  }
  isortit (dd, ind, numvertex); // ind stores indices in distance order

  //1 dump vertex

/**
   for (k=0; k<numvertex; ++k) {  
    kk = ind[k] ;
    node = &vlist[kk] ;
    freestring(&sss);  
    sss=strdup(node -> name) ;  
    substring(&sss, ":", "_") ;
    substring(&sss, "-", "_") ;
    fprintf(fff, "vertex %12s %9.0f\n", sss, node -> time) ;
   }
*/
  //2 dump label
  fprintf (fff, "digraph G { \n");
  fprintf (fff, "size = \"7.5,10\" ;\n");

  for (k = 0; k < numvertex; ++k) {
    kk = ind[k];
    node = &vlist[kk];
    if (node->label == NULL)
      continue;
    freestring (&sss);
    sss = strdup (node->name);
    substring (&sss, ":", "_");
    substring (&sss, "-", "_");
    fprintf (fff, "%12s ", sss);
    fprintf (fff, " [ label = \"%s\" ] ; ", node->label);

/**
    t = node -> popsize ;
    if ((totpop > 0) || (t>0)) {
     fprintf(fff, " %5d", t) ;
    }
*/
    fprintf (fff, "\n");
  }
// dump ledge, redge 
  vmax = 0.0;
  for (k = 0; k < numvertex; ++k) {
    kk = ind[k];
    node = &vlist[kk];
    edge = (EDGE *) node->eleft;
    if (edge != NULL) {
      xnode = (NODE *) node->left;
      vmax = MAX (vmax, edge->val);
      val = MAX (vmax * .0001, edge->val);
      pedge (fff, node, xnode, val, 0);
    }
    edge = (EDGE *) node->eright;
    if (edge != NULL) {
      xnode = (NODE *) node->right;
      vmax = MAX (vmax, edge->val);
      val = MAX (vmax * .0001, edge->val);
      pedge (fff, node, xnode, val, 0);
    }
  }
  for (k = 0; k < numvertex; ++k) {
    kk = ind[k];
    node = &vlist[kk];
    t = node->numwind;
    if (t == 0)
      continue;
    for (j = 0; j < t; ++j) {
      vind = node->windex[j];
      xnode = &vlist[vind];
      pedge (fff, xnode, node, node->wmix[j], 1);
    }
  }
  fprintf (fff, "} \n");
  fclose (fff);
  free (dd);
  free (ind);
  freestring (&sss);
}

void
pedge (FILE * fff, NODE * anode, NODE * bnode, double val, int mode)
{

  char s1[256], s2[256], s3[256];
  char *ss2 = NULL, *ss3 = NULL;
  int x;

  s1[0] = s2[0] = CNULL;

  if (mode == 1) {
    sprintf (s1, "%s", "style=dotted, ");
    x = nnint (100 * val);
    sprintf (s3, "%d%%", x);
  }

  else {
    x = nnint (1000 * val);
    sprintf (s3, "%d", x);
  }

  if (val > 0) {
    sprintf (s2, "label = \"%s\"", s3);
  }

  freestring (&ss2);
  ss2 = strdup (anode->name);
  substring (&ss2, ":", "_");
  substring (&ss2, "-", "_");
  fprintf (fff, "%s -> ", ss2);
  freestring (&ss3);
  ss3 = strdup (bnode->name);
  substring (&ss3, ":", "_");
  substring (&ss3, "-", "_");
  fprintf (fff, "%s ", ss3);
  fprintf (fff, "[ %s %s ] ; ", s1, s2);
  fprintf (fff, "\n");
  freestring (&ss2);
  freestring (&ss3);
}

void
dumpgraph (char *graphname)
{

  FILE *fff;
  NODE *node, *xnode, *xroot;
  EDGE *edge;
  int k, j, t, vind, kk;
  int *dd, *ind;
  char sform[10];

  if (hires)
    strcpy (sform, " %10.4f");
  else
    strcpy (sform, " %9.3f");
  if (graphname == NULL)
    fff = stdout;
  else
    openit (graphname, &fff, "w");

/** 
 we process vertices from the root 
*/
  xroot = root ();
  setdistances (xroot);
  ZALLOC (dd, numvertex, int);
  ZALLOC (ind, numvertex, int);
  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    dd[k] = node->distance;
  }
  isortit (dd, ind, numvertex); // ind stores indices in distance order

  //1 dump vertex
  for (k = 0; k < numvertex; ++k) {
    kk = ind[k];
    node = &vlist[kk];
    if (node->distance == HUGEDIS)
      fprintf (fff, "##");
    fprintf (fff, "vertex %12s %9.0f\n", node->name, node->time);
  }
  //2 dump label
  for (k = 0; k < numvertex; ++k) {
    kk = ind[k];
    node = &vlist[kk];
    if (node->isdead)
      continue;
    if (node->label == NULL)
      continue;
    fprintf (fff, "label %12s", node->name);
    fprintf (fff, " %20s", node->label);
    t = node->popsize;
    if ((totpop > 0) || (t > 0)) {
      fprintf (fff, " %5d", t);
    }
    fprintf (fff, "\n");
  }
// dump ledge, redge 
  for (k = 0; k < numvertex; ++k) {
    kk = ind[k];
    node = &vlist[kk];
    if (node->isdead)
      continue;
    edge = (EDGE *) node->eleft;
    if (edge != NULL) {
      xnode = (NODE *) node->left;
      fprintf (fff, "ledge %12s", edge->name);
      fprintf (fff, " %12s", node->name);
      fprintf (fff, " %12s", xnode->name);
      fprintf (fff, sform, edge->val);
      fprintf (fff, "\n");
    }
    edge = (EDGE *) node->eright;

    if (edge != NULL) {
      edge->val = MAX (edge->val, 0.0);
      xnode = (NODE *) node->right;
      fprintf (fff, "redge %12s", edge->name);
      fprintf (fff, " %12s", node->name);
      fprintf (fff, " %12s", xnode->name);
      fprintf (fff, sform, edge->val);
      fprintf (fff, "\n");
    }
  }
  for (k = 0; k < numvertex; ++k) {
    kk = ind[k];
    node = &vlist[kk];
    if (node->isdead)
      continue;
    t = node->numwind;
    if (t == 0)
      continue;
    fprintf (fff, "admix %12s  ", node->name);
    for (j = 0; j < t; ++j) {
      vind = node->windex[j];
      xnode = &vlist[vind];
      fprintf (fff, "%10s ", xnode->name);
    }
    for (j = 0; j < t; ++j) {
      fprintf (fff, sform, node->wmix[j]);
    }
    fprintf (fff, "\n");
  }
  if (graphname != NULL)
    fclose (fff);
  free (dd);
  free (ind);
}

double
exnames (char **pup, char **pdown, char **pename, int numedge)

/*up down edgename for edgenum */
{
  EDGE *edge;
  NODE *node;

  edge = &elist[numedge];
  *pename = strdup (edge->name);
  node = (NODE *) edge->up;
  *pup = strdup (node->name);
  node = (NODE *) edge->down;
  *pdown = strdup (node->name);
//printf("ename: %s  up:  %s\n", *pename, *pup) ;
  return (edge->val);

}

void
addnode (char *nodename, char *edgename, double breakval)
{

  NODE *node1, *node2, *newvert;
  EDGE *edge, *newedge;
  int t;
  char ss[MAXSTR];


  t = edgenum (edgename);
  if (t < 0)
    fatalx ("edge: %s not found\n", edgename);
  edge = &elist[t];
  node1 = (NODE *) edge->up;
  node2 = (NODE *) edge->down;
  newedge = &elist[numedge];
  newvert = &vlist[numvertex];
  ++numedge;
  ++numvertex;
  newvert->name = strdup (nodename);

  strcpy (ss, edge->name);
  strcat (ss, ":2");
  newedge->name = strdup (ss);

  strcpy (ss, edge->name);
  strcat (ss, ":1");
  freestring (&edge->name);
  edge->name = strdup (ss);

  newedge->val = edge->val * (1.0 - breakval);
  edge->val *= breakval;

  if (edge->isleft) {
    node1->left = (struct NODE *) newvert;
  }
  else {
    node1->right = (struct NODE *) newvert;
  }
  newvert->eleft = (struct EDGE *) newedge;
  newvert->left = (struct NODE *) node2;
  newedge->isleft = YES;
  newedge->down = (struct NODE *) node2;
  newedge->up = (struct NODE *) newvert;
  edge->down = (struct NODE *) newvert;
  newvert->parent = (struct NODE *) node1;
  node2->parent = (struct NODE *) newvert;

  newvert->eparent = (struct EDGE *) edge;
  node2->eparent = (struct EDGE *) newedge;

}

void
addmixedge2 (char **nodes, char **edgenames, char *newedgename,
             char *label, double *breakval, double eval, double *wmix,
             double *eeval, int nadmix)
// label can be NULL
// xtra edge introduced... 
{
  NODE *node, *newvert;
  EDGE *edge, *newedge;
  char *nodename, *n2name, **nnodes;
  int t, k;
  char s1[MAXSTR], s2[MAXSTR];

  ZALLOC (nnodes, nadmix, char *);
  for (k = 0; k < nadmix; ++k) {
    addnode (nodes[k], edgenames[k], breakval[k]);
    sprintf (s1, "%s:3", nodes[k]);
    nnodes[k] = strdup (s1);
    sprintf (s2, "%s:3", edgenames[k]);
    addedge (s2, nodes[k], nnodes[k]);
    t = edgenum (s2);
    newedge = &elist[t];
    newedge->val = eeval[k];
  }

  nodename = nodes[nadmix];
  n2name = nodes[nadmix + 1];

  addadmix (nodename, nnodes, wmix, nadmix);
  addedge (newedgename, nodename, n2name);

  t = edgenum (newedgename);
  newedge = &elist[t];
  newedge->val = eval;

  if (label != NULL) {
    t = vertexnum (n2name);
    newvert = &vlist[t];
    newvert->label = strdup (label);
    eglist[numpops] = strdup (label);
    egnum[numpops] = t;
    ++numpops;
  }
  freeup (nnodes, nadmix);
}

void
addedgenode (char *nodename, char *n2name, char *edgename, char *newedgename,
             char *label, double breakval, double eval)
// label can be NULL
{
  NODE *node, *newvert;
  EDGE *edge, *newedge;
  int t;


  addnode (nodename, edgename, breakval);
  addedge (newedgename, nodename, n2name);
  t = edgenum (newedgename);
  newedge = &elist[t];
  newedge->val = eval;

  if (label != NULL) {
    t = vertexnum (n2name);
    newvert = &vlist[t];
    newvert->label = strdup (label);
    eglist[numpops] = strdup (label);
    egnum[numpops] = t;
    ++numpops;
  }
}

void
addmixedge (char **nodes, char **edgenames, char *newedgename,
            char *label, double *breakval, double eval, double *wmix,
            int nadmix)
// label can be NULL
{
  NODE *node, *newvert;
  EDGE *edge, *newedge;
  char *nodename, *n2name;
  int t, k;
  char s1[MAXSTR], s2[MAXSTR];

  for (k = 0; k < nadmix; ++k) {
    addnode (nodes[k], edgenames[k], breakval[k]);
  }

  nodename = nodes[nadmix];
  n2name = nodes[nadmix + 1];

  addadmix (nodename, nodes, wmix, nadmix);
  addedge (newedgename, nodename, n2name);
  t = edgenum (newedgename);
  newedge = &elist[t];
  newedge->val = eval;

  if (label != NULL) {
    t = vertexnum (n2name);
    newvert = &vlist[t];
    newvert->label = strdup (label);
    eglist[numpops] = strdup (label);
    egnum[numpops] = t;
    ++numpops;
  }
}

int
vertexnum (char *vertname)
{
  NODE *node;
  int t, k;

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    t = strcmp (vertname, node->name);
    if (t == 0)
      return k;
  }
  return -1;
}

int
edgenum (char *edgename)
{
  EDGE *edge;
  int t, k;

  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    t = strcmp (edgename, edge->name);
    if (t == 0)
      return k;
  }
  return -1;

}

int
setxx (NODE * node, NODE ** xx)
{
  int t = 0;
  NODE *nn;

  nn = (NODE *) node->parent;
  if (nn != NULL) {
    xx[t] = nn;
    ++t;
  }
  nn = (NODE *) node->left;
  if (nn != NULL) {
    xx[t] = nn;
    ++t;
  }
  nn = (NODE *) node->right;
  if (nn != NULL) {
    xx[t] = nn;
    ++t;
  }

  return t;

}

void
reroot (char *nodename)
{

  NODE *node, *oldroot, *tnode, *nn[4], *ttnode, *xx[3];
  EDGE *edge, *newedge, ee[3];
  int nee, nset, t, dd[4], nw, k, j, unset, dmin, jmin;
  int numiter;

  node = vert (nodename);
  if (node == NULL)
    fatalx ("node: %s not found\n", nodename);
  oldroot = root ();
  printf ("oldroot: %s\n", oldroot->name);
  nee = setxx (node, xx);       // neighbours
  if (nee == 3)
    fatalx ("3 edges from new root: %s\n", nodename);

/* set distances */
  setdistances (node);
  oldroot->isroot = NO;
  node->isroot = NO;
  for (k = 0; k < numvertex; ++k) {
    tnode = &vlist[k];
    fixnode (tnode);
  }
  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    fixedge (edge);
  }
  node->isroot = YES;
}

void
setdistances (NODE * node)
{
  NODE *oldroot, *tnode, *nn[4], *ttnode, *xx[3], *unode;
  EDGE *edge, *newedge, ee[3];
  int nee, nset, t, dd[4], nw, k, j, unset, dmin, jmin;
  int numiter;

  cleardis ();
  oldroot = root ();
  setdis (node, 0);
  if (oldroot->distance > 1000)
    fatalx ("%s must be in base tree\n", node->name);
  numiter = 0;
  for (;;) {
    nset = 0;
    ++numiter;
    if (numiter > 20)
      fatalx ("looping...\n");
    unset = NO;                 // vertex unset ??  
    for (k = 0; k < numvertex; ++k) {
      tnode = &vlist[k];
      if (tnode->isdead)
        continue;
      nw = tnode->numwind;
      if (nw == 0)
        continue;
      if (tnode->distance < HUGEDIS)
        continue;
      for (j = 0; j < nw; ++j) {
        ttnode = nn[j] = &vlist[tnode->windex[j]];
        dd[j] = ttnode->distance;
      }
      ivlmaxmin (dd, nw, NULL, &jmin);
      dmin = dd[jmin];
      if (dmin > 999999) {
        unset = YES;
        unode = nn[jmin];
        continue;               // not set
      }
      setdis (tnode, dmin + 1000);      // next tree
      ++nset;
    }
    if (nset == 0)
      break;
  }
  if (unset)
    fatalx ("node not linked to root... %s\n", unode->name);
}

void
fixedge (EDGE * edge)
{
  NODE *n1, *n2;

  n1 = (NODE *) edge->up;
  n2 = (NODE *) edge->down;

  if (n1->distance < n2->distance)
    return;
  edge->up = (struct NODE *) n2;
  edge->down = (struct NODE *) n1;
}

void
fixnode (NODE * node)
// fix up edges 
{
  int qleft = NO, j, nee, basedis, ddis;
  NODE *xx[3], *tnode;

  if (node->numwind > 0)
    return;
  nee = setxx (node, xx);       // neighbours
  node->parent = node->left = node->right = NULL;
  basedis = node->distance;
  for (j = 0; j < nee; ++j) {
    tnode = xx[j];
    ddis = tnode->distance;
    if (ddis < basedis) {
      node->parent = (struct NODE *) tnode;
      continue;
    }
    if (qleft)
      node->right = (struct NODE *) tnode;
    if (!qleft)
      node->left = (struct NODE *) tnode;
    qleft = YES;
  }
  node->eparent = (struct EDGE *) xedge (node, (NODE *) node->parent);
  if (node->eparent != NULL) {
    tnode = (NODE *) node->parent;
//  printf("qq2 %s\n", tnode -> name) ; 
  }
  node->eleft = (struct EDGE *) xedge (node, (NODE *) node->left);
  node->eright = (struct EDGE *) xedge (node, (NODE *) node->right);
}

EDGE *
xedge (NODE * n1, NODE * n2)
// find edge linking two vertices
{
  EDGE *edge;
  int k;
  NODE *x1, *x2;

  if (n1 == NULL)
    return NULL;
  if (n2 == NULL)
    return NULL;
  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    x1 = (NODE *) edge->up;
    x2 = (NODE *) edge->down;
    if ((n1 == x1) && (n2 == x2))
      return edge;
    if ((n1 == x2) && (n2 == x1))
      return edge;
  }
  return NULL;
}

void
setdis (NODE * node, int dis)
{
  NODE *xx[3];
  int nee, j;

  if (node->distance < dis)
    return;
  node->distance = MIN (node->distance, dis);
  nee = setxx (node, xx);       // neighbours
  for (j = 0; j < nee; ++j) {
    setdis (xx[j], dis + 1);
  }
}

void
cleardis ()
{
  NODE *node;
  int t, k;

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    node->distance = HUGEDIS;
  }

}

NODE *
vert (char *nodename)
{
  NODE *node;
  int t, k;

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    t = strcmp (nodename, node->name);
    if (t == 0)
      return node;
  }
  return NULL;

}

NODE *
root ()
{
  NODE *node;
  int t, k;

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    if (node->numwind > 0)
      continue;
    if (node->eparent != NULL)
      continue;
    return node;
  }
  fatalx ("can't find root\n");
}

void
settime ()
// fill in times on nodes
{
  double konst = 10000;         // nominal scaling for drift
  double *eq, *rhs, *pp, *ans;
  int nneq, neq, nv, t, v1, v2, j, x, k;
  int ww[100];
  NODE *xroot, *n1, *n2, *node;
  EDGE *edge;
  int *vnum;

  nneq = 2 * (numedge + numvertex) + 10000;
  nv = numvertex;

/**
  xroot = root() ;
  reroot(xroot) ;
*/

  ZALLOC (eq, nneq * nv, double);
  ZALLOC (rhs, nneq, double);
  ZALLOC (vnum, nv, int);
  ZALLOC (ans, nv, double);
  idperm (vnum, nv);
  setvnum (vnum, egnum, numpops);
  for (j = 0; j < numvertex; ++j) {
    node = &vlist[j];
    t = node->numwind;
    if (t == 0)
      continue;
    setvnum (vnum, node->windex, t);
  }
  printimat (vnum, 1, nv);

  pp = eq;
  for (k = 0; k < nv; ++k) {
    n1 = &vlist[k];
    pp[k] = 1.0e-6;
    if (n1->isroot)
      pp[k] = 1.0e6;
    pp += nv;
  }
  neq = nv;
  for (j = 0; j < numedge; ++j) {
    edge = &elist[j];
    n1 = (NODE *) edge->up;
    n2 = (NODE *) edge->down;
    rhs[neq] = edge->val * konst;
    v1 = n1->gnode;
    v2 = n2->gnode;
    v1 = vnum[v1];
    v2 = vnum[v2];
    printf ("zz %s %d %d\n", n1->name, n1->gnode, v1);
    printf ("zz %s %d %d\n", n2->name, n2->gnode, v2);
    pp[v2] = 1.0;
    pp[v1] = -1.0;
    pp += nv;
    ++neq;
  }
  printmat (eq, neq, nv);
  printnl ();
  printnl ();
  printf ("rhs:\n");
  printmat (rhs, 1, neq);
  printnl ();
  regressit (ans, eq, rhs, neq, nv);
  printf ("ans:\n");
  printmat (ans, 1, nv);
  printnl ();
  for (k = 0; k < nv; ++k) {
    v1 = vnum[k];
    n1 = &vlist[k];
    n1->time = ans[v1];
  }
  free (eq);
  free (rhs);
  free (vnum);
}

void
setvnum (int *vnum, int *list, int n)
{

  int a, k, x;

  x = 999999999;
  for (a = 0; a < n; a++) {
    k = list[a];
    x = MIN (x, vnum[k]);
  }
  for (a = 0; a < n; a++) {
    k = list[a];
    vnum[k] = x;
  }
}


int
qreadit (char *cname)
// main input routine for graph (new format)
{
#define MAXSTR 128
#define MAXFF  10
  FILE *fff;
  char line[MAXSTR + 1], c;
  char *spt[MAXFF], *sx, *s1, *s2;
  int nsplit, n = 0;
  int kret, t, v1, v2, w, k, nv;
  int okline;
  NODE *node, *tnode, *node2, *node1;
  EDGE *edge;
  double ww[4];
  char *ssx[4];
  int nt;
  double val;


  openit (cname, &fff, "r");
  line[MAXSTR] = '\0';
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    okline = NO;

    kret = strcmp (sx, "root");
    if (kret == 0) {
      trivgraph (spt[1]);
      okline = YES;
    }
    kret = strcmp (sx, "edge");
    if (kret == 0) {
      ckline (line, nsplit, 4);
      addedge (spt[1], spt[2], spt[3]);
      okline = YES;
    }
    kret = strcmp (sx, "label");
    if (kret == 0) {
      ckline (line, nsplit, 3);
      addlabel (spt[1], spt[2]);
      okline = YES;
    }
    kret = strcmp (sx, "admix");
    if (kret == 0) {
      okline = YES;
      getadwts (spt, ssx, ww, nsplit, &nt);
      addadmix (spt[1], ssx, ww, nt);
      okline = YES;
    }

    kret = strcmp (sx, "lock");
    if (kret == 0) {
      okline = YES;
      nv = getnumvertex ();
      k = vindex (spt[1], vlist, nv);
      if (k < 0)
        fatalx ("bad line %d", line);
      tnode = &vlist[k];
      tnode->isfixed = YES;
    }
    kret = strcmp (sx, "edgelock");
    if (kret == 0) {
      okline = YES;
      k = edgenum (spt[1]);
      if (k < 0)
        fatalx ("bad line %d", line);
      edge = &elist[k];
      val = atof (spt[2]);
      edge->val = atof (spt[2]);
      edge->isfixed = YES;
    }


    if (okline == NO) {
      printf ("*** warning bad line:\n %s\n", line);
    }
    freeup (spt, nsplit);
    continue;
  }


  fclose (fff);
  return numvertex;

}

void
ckline (char *line, int ns, int n)
{
  if (ns < n)
    fatalx ("bad line %s\n", line);
}

void
superalloc (int **xv, int ***xe, int ***adv, int **aedge)
{

  int t;

  ZALLOC (*xv, numvertex, int);
  ZALLOC (*aedge, numadmix, int);

  *xe = initarray_2Dint (numedge, 3, -1);
  t = MAX (numadmix, 1);
  *adv = initarray_2Dint (t, MAXW, -1);
}

void
supersetup (int *xvlist, int *nxvlist, int **xelist, int *nxelist,
            int **admixv, int *admixedge, int *nxalist)
{
  NODE *xroot, *n1, *n2, *node;
  EDGE *edge;
  int *vnum;
  int nxv, nxe, nxa;
  int k, t;


  fflush (stdout);
  nxa = nxv = nxe = 0;

  for (k = 0; k < numedge; ++k) {
    edge = &elist[k];
    n1 = (NODE *) edge->up;
    n2 = (NODE *) edge->down;
    xelist[nxe][0] = n1->gnode;
    xelist[nxe][1] = n2->gnode;
    xelist[nxe][2] = edge->edgenum;
    ++nxe;
  }
  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    if (node->isadmix)
      continue;
    node->supernum = nxv;
    xvlist[nxv] = k;
    ++nxv;
  }
  for (k = 0; k < numadmix; ++k) {
    node = adlist[k];
    ivclear (admixv[nxa], -1, MAXW);
    t = node->numwind;
    if (t <= 0)
      fatalx ("bad logic\n");
    copyiarr (node->windex, admixv[nxa], t);

    if (node->left == NULL)
      fatalx ("bad admix node %s\n", node->name);
    if (node->right != NULL)
      fatalx ("bad admix node (redge) %s\n", node->name);
    edge = (EDGE *) node->eleft;
    admixedge[nxa] = edge->edgenum;
    printf ("zz %s %d\n", edge->name, edge->edgenum);
    fflush (stdout);
    ++nxa;
  }

  *nxvlist = nxv;
  *nxelist = nxe;
  *nxalist = nxa;
}

void
superest (double *xmean, double *xvar, double *svar, double *yobs,
          int nxvlist, int *eelist)
{
  double tiny = 1.0e-12;
  int x, a, t;
  int nv2 = nxvlist * nxvlist;
  double *v1, *v2, *rhs, *ww;
  int *xtype;
  NODE *node;

  ZALLOC (v1, nv2, double);
  ZALLOC (v2, nv2, double);
  ZALLOC (rhs, nxvlist, double);
  ZALLOC (ww, nxvlist, double);
  ZALLOC (xtype, nxvlist, int);

  copyarr (svar, v1, nv2);
  for (x = 0; x < nxvlist; ++x) {
    v1[x * nxvlist + x] += tiny;
  }
  pdinv (v2, v1, nxvlist);
  for (x = 0; x < numpops; ++x) {
    t = eelist[x];
    rhs[t] = ww[t] = yobs[x];
    xtype[t] = 1;
  }
  for (a = 0; a < nxvlist; ++a) {
    if (xtype[a] == 1)
      continue;
    rhs[a] = -vdot (v2 + a * nxvlist, ww, nxvlist);
  }
  for (x = 0; x < numpops; ++x) {
    t = eelist[x];
    for (a = 0; a < nxvlist; ++a) {
      v2[t * nxvlist + a] = v2[a * nxvlist + t] = 0.0;
    }
    v2[t * nxvlist + t] = 1.0;
  }
  pdinv (v1, v2, nxvlist);
// mulmat(xmean, v1, rhs, nxvlist, nxvlist, 1) ;
  solvit (v2, rhs, nxvlist, xmean);
  copyarr (v1, xvar, nv2);

  for (x = 0; x < numpops; ++x) {
    t = eelist[x];
    for (a = 0; a < nxvlist; ++a) {
      xvar[t * nxvlist + a] = xvar[a * nxvlist + t] = 0.0;
    }
  }

  free (v1);
  free (v2);
  free (rhs);
  free (ww);
  free (xtype);

}

void
supergetvnames (char **vnames, int *xvlist, int nxvlist)
{
  int a, x;
  NODE *node;

  for (a = 0; a < nxvlist; ++a) {
    x = xvlist[a];
    node = &vlist[x];
    vnames[a] = node->name;
  }

}

void
supergetvar (double *svar,
             int *xvlist, int nxvlist, int **xelist, int nxelist,
             int **admixv, int *admixedge, int nxalist)
{
  double *pwts;
  double *ee, *ww;
  char **vnames;
  int t, a, b, k, x1, x2, nrows, ncols;
  double y;
  NODE *node;
  EDGE *edge;

  ZALLOC (pwts, numvertex * numedge, double);
  ZALLOC (vnames, numvertex, char *);
  ZALLOC (ee, numedge, double);
  ZALLOC (ww, numedge, double);

  getpnwts (pwts, &nrows, &ncols);
  vzero (svar, nxvlist * nxvlist);

  supergetvnames (vnames, xvlist, nxvlist);
  for (a = 0; a < nxelist; ++a) {
    k = xelist[a][2];
    edge = &elist[k];
    t = edge->edgenum;
    ee[t] = edge->val;
  }

  for (a = 0; a < nxvlist; ++a) {
    x1 = xvlist[a];
    for (b = 0; b < nxvlist; ++b) {
      x2 = xvlist[b];

      vvt (ww, pwts + x1 * numedge, pwts + x2 * numedge, numedge);
      vvt (ww, ww, ee, numedge);
      y = asum (ww, numedge);

      svar[a * nxvlist + b] = y;
      svar[b * nxvlist + a] = y;
    }
  }

/**
  printf("supergetvar: \n") ;
  printmatz(svar, vnames, nxvlist) ;
  printnl() ;
*/

  free (ee);
  free (ww);
  free (pwts);
  free (vnames);

}

void
superreest (double *s2,
            int *xvlist, int nxvlist, int **xelist, int nxelist, int **admixv,
            int *admixedge, int nxalist)
{
#define MAXP (MAXW+1)
  int i, xn1, xn2, z1, z2, eenum, k, t, tp, kk, a, b;
  double vv[MAXP], ww[MAXP], xx[MAXP * MAXP], xxx[MAXW * MAXW];
  double y;
  NODE *n1, *n2, *nn;
  EDGE *edge;
  int cc[MAXP];
  static int ncall = 0;

  ++ncall;
  for (i = 0; i < nxelist; ++i) {
    xn1 = xelist[i][0];
    xn2 = xelist[i][1];
    eenum = xelist[i][2];
    n1 = &vlist[xn1];
    n2 = &vlist[xn2];
    edge = &elist[eenum];
    if (n1->isadmix)
      continue;

    z1 = n1->supernum;
    z2 = n2->supernum;

    cc[0] = z1;
    cc[1] = z2;
    ww[0] = 1.0;
    ww[1] = -1.0;
    t = 1;
    rcsquish (xx, s2, cc, nxvlist, t + 1);
    mulmat (vv, xx, ww, t + 1, t + 1, 1);
    y = vdot (vv, ww, t + 1);
    edge->val = y;

/**
   y  = s2[z1*nxvlist+z1] ;  
   y += s2[z2*nxvlist+z2] ;  
   y -= 2.0*s2[z1*nxvlist+z2] ;
   edge -> val = y ;
*/

  }

  for (i = 0; i < nxelist; ++i) {
    break;
    xn1 = xelist[i][0];
    xn2 = xelist[i][1];
    eenum = xelist[i][2];
    n1 = &vlist[xn1];
    n2 = &vlist[xn2];
    edge = &elist[eenum];
    if (n1->isadmix == NO)
      continue;
    t = n1->numwind;
    for (k = 0; k < t; ++k) {
      kk = n1->windex[k];
      nn = &vlist[kk];
      cc[k] = nn->supernum;
      ww[k] = n1->wmix[k];
    }
    bal1 (ww, t);
    cc[t] = n2->supernum;
    ww[t] = -1.0;
    rcsquish (xx, s2, cc, nxvlist, t + 1);
    mulmat (vv, xx, ww, t + 1, t + 1, 1);
    y = vdot (vv, ww, t + 1);
    edge->val = y;
  }
// weight estimation
  for (i = 0; i < nxelist; ++i) {
    xn1 = xelist[i][0];
    xn2 = xelist[i][1];
    eenum = xelist[i][2];
    n1 = &vlist[xn1];
    n2 = &vlist[xn2];
    edge = &elist[eenum];
    if (n1->isadmix == NO)
      continue;
    t = n1->numwind;
    for (k = 0; k < t; ++k) {
      kk = n1->windex[k];
      nn = &vlist[kk];
      cc[k] = nn->supernum;
    }

    cc[t] = n2->supernum;
    tp = t + 1;
    rcsquish (xx, s2, cc, nxvlist, tp);
    copyarr (n1->wmix, ww, t);

    ww[t] = -1.0;
    if (ncall < 10) {
      printf ("zz %s %d ", n1->name, n1->isfixed);
      printmat (n1->wmix, 1, t);
    }

    if (n1->isfixed == YES) {
      mulmat (vv, xx, ww, t + 1, t + 1, 1);
      y = vdot (vv, ww, t + 1);
      edge->val = y;
      continue;
    }
    for (a = 0; a < t; a++) {
      for (b = 0; b < t; b++) {
        xxx[a * t + b] = xx[a * tp + b] + xx[t * tp + t];
        xxx[a * t + b] -= xx[a * tp + t];
        xxx[a * t + b] -= xx[b * tp + t];
      }
    }
    y = qwmax (xxx, ww, t);

/**
   printf("zzqw\n") ;
   printmat(xxx, t, t) ;
   printnl() ;
   printmat(ww, 1,  t) ;
   printf("zzqwans %12.6f %12.6f\n", edge -> val, y) ;
*/
    copyarr (ww, n1->wmix, t);
    edge->val = y;
    if (ncall < 10) {
      printf ("zz %s %d ", n1->name, n1->isfixed);
      printmat (n1->wmix, 1, t);
    }
  }

  return;

}

void
supergetvals (double **admixw, double *elen,
              int *xvlist, int nxvlist, int **xelist, int nxelist,
              int **admixv, int *admixedge, int nxalist)
{
  double *ewts, **vmix;
  int k, t, eenum;
  EDGE *edge;
  NODE *node;

  ZALLOC (ewts, numedge, double);
  vmix = initarray_2Ddouble (MAXG, MAXW, 0.0);

  getewts (ewts);
  for (k = 0; k < nxelist; ++k) {
    t = xelist[k][2];
    elen[k] = ewts[t];
  }

  for (k = 0; k < nxalist; ++k) {
    eenum = admixedge[k];
    edge = &elist[eenum];
    node = (NODE *) edge->up;
    t = node->numadmix;
    if (t < 0)
      fatalx ("logic bug: %s %s\n", node->name, edge->name);
    copyarr (node->wmix, admixw[k], t);
  }
}

void
superputvals (double **admixw, double *elen,
              int *xvlist, int nxvlist, int **xelist, int nxelist,
              int **admixv, int *admixedge, int nxalist)
{
  double *ewts, **vmix;
  int k, t, eenum;
  EDGE *edge;
  NODE *node;

  ZALLOC (ewts, numedge, double);
  vmix = initarray_2Ddouble (MAXG, MAXW, 0.0);

  for (k = 0; k < nxelist; ++k) {
    t = xelist[k][2];
    ewts[t] = elen[k];
  }
  putewts (ewts);

  for (k = 0; k < nxalist; ++k) {
    eenum = admixedge[k];
    edge = &elist[eenum];
    node = (NODE *) edge->up;
    t = node->numadmix;
    if (t < 0)
      fatalx ("logic bug: %s %s\n", node->name, edge->name);
    copyarr (admixw[k], vmix[t], MAXW);
    if (node->isfixed) {
      copyarr (node->wmix, vmix[t], MAXW);
    }
  }
  putgmix (vmix);


  free (ewts);
  free2D (&vmix, MAXG);

}

void
supereglist (int *eelist)
{
  int k, t;
  NODE *node;

  ivzero (eelist, numvertex);
  for (t = 0; t < numpops; ++t) {
    k = egnum[t];
    node = &vlist[k];
    eelist[t] = node->supernum;
  }
}

void
setadmfix (char *fixname)
{
  FILE *fff;
  char line[MAXSTR + 1], c;
  char *spt[MAXFF], *sx, *s1, *s2;
  int nsplit, n = 0;
  int kret, t, a, w, k;
  int okline;
  NODE *node, *tnode, *node2, *node1;
  EDGE *edge;
  double ww[4];
  char *ssx[4];
  int nv;
  double val;


  openit (fixname, &fff, "r");
  printf ("setadmfix called\n");
  fflush (stdout);
  line[MAXSTR] = '\0';
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    nv = getnumvertex ();
    k = vindex (sx, vlist, nv);

    if (k < 0) {

      for (a = 0; a < nv; a++) {
        node = &vlist[a];
        printf ("zzbad %3d %s %s\n", a, sx, node->name);
      }

      fatalx ("bad fix line: %s\n", line);
    }

    node = &vlist[k];
    t = node->numwind;
    for (a = 0; a < t; ++a) {
      node->wmix[a] = atof (spt[a + 1]);
    }
    bal1 (node->wmix, t);
    node->isfixed = YES;
    printf ("node %s fixed: ", node->name);
    printmat (node->wmix, 1, t);
    fflush (stdout);
  }
  fclose (fff);
  printf ("setfix done\n");
  fflush (stdout);
}

int
findlabel (char *label)
{
  int k, t = -1;
  NODE *node;
  char *lab;

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];
    lab = node->label;
    if (lab == NULL)
      continue;
    if (strcmp (lab, label) == 0)
      return k;
  }
  return -1;

}

void
getpops (char **eelist, int *npops)
// eelist preallocated
{
  *npops = numpops;
  copystrings (eglist, eelist, numpops);
  return;
}

void
copystringsd (char **eglist, char **neweglist, int numeg, int xdel)
{
  int i, j;

  if (xdel < 0)
    fatalx ("yukk\n");
  if (xdel >= numeg)
    fatalx ("yukk\n");
  j = 0;
  for (i = 0; i < numeg; ++i) {
    if (i == xdel)
      continue;
    neweglist[j] = strdup (eglist[i]);
    ++j;
  }
}

int
hashgraph ()
{
  int *hashvals;
  NODE *node, *node2;
  int k, t, rootnum;

  ZALLOC (hashvals, numvertex, int);

  for (k = 0; k < numvertex; ++k) {
    node = &vlist[k];

    if (node->label != NULL) {
      t = stringhash (node->label);
      hashvals[k] += t;
    }
    if (node->isleaf) {
      hashg (k, hashvals);
    }
  }
  node = root ();
  rootnum = node->gnode;

  if (rootnum < 0)
    fatalx ("no root\n");

  t = hashvals[rootnum];

  free (hashvals);
  return t;
}

void
hashg (int knum, int *hashvals)
{
  NODE *node, *nodep;
  int a, t, kk;

  node = &vlist[knum];

  t = node->numwind;
  for (a = 0; a < t; a++) {
    kk = node->windex[a];
    hashvals[kk] += (hashvals[knum] + 1001);
    hashg (kk, hashvals);
  }
  nodep = (NODE *) node->parent;
  if (nodep != NULL) {
    kk = nodep->gnode;
    hashvals[kk] += (hashvals[knum] + 1);
    hashg (kk, hashvals);
  }
}
