HOMEL=$(PWD)
TOP=../bin
BIN=$(HOMEL)/../bin

override LDLIBS += -lgsl -lopenblas -lm -lnick
# Some Linux distributions require separate lapacke library
# override LDLIBS += -llapacke
override LDFLAGS += -g  -p -pg -L./nicksrc
override CFLAGS += -c -g -p -pg -Wimplicit -I./ -I./nicksrc
# Mac additions using homebrew installations
#override CFLAGS += -I/usr/local/opt/openblas/include -I/usr/local/opt/gsl/include
#override LDFLAGS += -L/usr/local/opt/openblas/lib -L/usr/local/opt/gsl/lib
# Harvard Medical School O2 cluster additions
ifdef SLURM_CONF
override CFLAGS += -I/n/app/openblas/0.2.19/include -I/n/app/gsl/2.3/include
override LDFLAGS += -L/n/app/openblas/0.2.19/lib -L/n/app/gsl/2.3/lib/
TOP=../bin
BIN=$(HOMEL)/../bin
endif

ND = nicksrc
NLIB = $(ND)/libnick.a

PROGS= qp3Pop qpDstat qpF4ratio qpAdm qpWave qp4diff dowtjack qpBound qpGraph qpreroot qpff3base qpDpart qpfstats qpfmv  qpmix 
PROGS2 = rexpfit.r  wtjack.pl 
PROGS3 = convertf mergeit snpunion simpjack2 grabpars  easystats easycheck easylite  multimerge geno_single transpose merge_transpose
PROGS4 = multimerge merge_transpose transpose geno_single nickhash
PERLP = jackdiff mkpretty numlines getresult xtractcol 

CC = cc

all: $(NLIB) $(PROGS) $(PROGS3)  $(PROGS4) 


$(NLIB):
	$(MAKE) -C $(ND)

## an ubuntu user found he needed: make all LDLIBS="-llapack" 

statsubs.o:     tables nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -o statsubs.o nicksrc/statsubs.c


qpDstat: $(NLIB) qpDstat.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qpdslow: $(NLIB) qpdslow.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qpDpart: $(NLIB) qpDpart.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qp3Pop: $(NLIB) qp3Pop.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qp4diff: $(NLIB) qp4diff.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qpF4ratio: $(NLIB) qpF4ratio.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qpF4ratio: $(NLIB) qpF4ratio.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qpf4diff: $(NLIB) f4diff.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

easystats: $(NLIB) easystats.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o
easycheck: $(NLIB) easycheck.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o
easylite: $(NLIB) easylite.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

multimerge: multimerge.o mcio.o admutils.o egsubs.o geno.o $(NLIB)
nickhash:  nickhash.o admutils.o $(NLIB)
merge_transpose: merge_transpose.o geno.o
transpose: transpose.o geno.o
geno_single: geno_single.o geno.o

rolloff: $(NLIB) rolloff.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

rolloffp: $(NLIB) rolloffp.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qpBound: $(NLIB) qpBound.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o geno.o

dowtjack: $(NLIB) dowtjack.o

convertf: $(NLIB) convertf.o mcio.o mapfile.o admutils.o egsubs.o h2d.o exclude.o  geno.o

snpunion: $(NLIB) snpunion.o mcio.o mapfile.o admutils.o egsubs.o h2d.o exclude.o  geno.o

mergeit: $(NLIB) mergeit.o mcio.o mapfile.o admutils.o egsubs.o h2d.o exclude.o  geno.o

grabpars: $(NLIB) grabpars.o

qpff3base: $(NLIB) qpff3base.o qpgsubs.o mapfile.o qpsubs.o mcio.o admutils.o egsubs.o regsubs.o ldsubs.o geno.o

qpfstats: $(NLIB) qpfstats.o qpgsubs.o mapfile.o qpsubs.o mcio.o admutils.o egsubs.o regsubs.o ldsubs.o geno.o

qqpfstats: $(NLIB) qqpfstats.o qpgsubs.o qpsubs.o mcio.o mapfile.o admutils.o egsubs.o regsubs.o ldsubs.o

qpfmv: $(NLIB) qpfmv.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o  eigsubs.o eigx.o geno.o

qpfmvmix: $(NLIB) qpfmvmix.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o  eigsubs.o eigx..  geno.o

qqpfmv: $(NLIB) qqpfmv.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o  eigsubs.o eigx.o

qpAdm: $(NLIB) qpAdm.o qpsubs.o mcio.o ldsubs.o mapfile.o admutils.o  egsubs.o regsubs.o f4rank.o  eigsubs.o eigx.o geno.o
qpfit: $(NLIB) qpfit.o qpsubs.o mcio.o ldsubs.o mapfile.o admutils.o  egsubs.o regsubs.o f4rank.o  eigsubs.o eigx.o geno.o

qqpAdm: $(NLIB) qqpAdm.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o f4rank.o  eigsubs.o eigx.o
qqAdm: $(NLIB) qqAdm.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o f4rank.o  eigsubs.o eigx.o

qpWave: $(NLIB) qpWave.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o f4rank.o eigsubs.o eigx.o geno.o

qpmix:  $(NLIB)  qpmix.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o eigsubs.o eigx.o regsubs.o qpfitr1.o gslrank1.o  geno.o
qptestr1:  $(NLIB)  qptestr1.o  gslrank1.o 

oldqpmix:  $(NLIB)  oldqpmix.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o
qqpmix:  $(NLIB)  qqpmix.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o
mkweights:  $(NLIB)  mkweights.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o
 
w2p:  $(NLIB)  w2p.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o  egsubs.o regsubs.o

simpjack2: $(NLIB) simpjack2.o

install:	all dirs tables
	mkdir -p $(TOP)
	cp $(PROGS) $(TOP)
	cp $(PROGS2) $(TOP)
	cp $(PROGS3) $(TOP)
	cp $(PROGS3) $(TOP)
	cp ../perlsrc/* $(TOP) 
	cp ../pysrc/* $(TOP) 
	cp fxtract $(BIN)
	cp -r script $(BIN) 

uninstall:
	rm -f tables
	rm -f $(TOP)/*

tables:    
	echo "tables made"  > tables
	cp twtable  $(HOMEL)/admixtables
	
dirs:	
	mkdir -p  $(HOMEL)/admixtables

clean: 
	rm -f *.o 
	rm -f *junk*
	rm -f core
	$(MAKE) -C $(ND) clean
	rm -f $(PROGS) $(PROGS3) $(PROGS4) 

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/admixtables 
	rm -rf $(BIN)/script 

## gsl-config may help here


qpGraph: $(NLIB) qpGraph.o gslqp.o qpgsubs.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

qqpGraph: $(NLIB) qqpGraph.o gslqp.o qpgsubs.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpreroot: qpreroot.o qpgsubs.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o geno.o

diffmean: diffmean.o 

gcount: gcount.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o  geno.o 

kimf: kimf.o gslkim.o qpgsubs.o qpsubs.o mcio.o mapfile.o ldsubs.o admutils.o egsubs.o regsubs.o

