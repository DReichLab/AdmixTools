HOMEL=$(PWD)
TOP=../bin
BIN=$(HOMEL)/../bin

override LDLIBS += -lgsl -lopenblas -lm -lnick
# Some Linux distributions require separate lapacke library
# override LDLIBS += -llapacke
override LDFLAGS += -L./nicksrc
override CFLAGS += -c -g -p -Wimplicit -I./ -I./nicksrc
# Mac additions using homebrew installations
#override CFLAGS += -I/usr/local/opt/openblas/include -I/usr/local/opt/gsl/include
#override LDFLAGS += -L/usr/local/opt/openblas/lib -L/usr/local/opt/gsl/lib
# Harvard Medical School O2 cluster additions
ifdef SLURM_CONF
override CFLAGS += -I/n/app/openblas/0.2.19/include -I/n/app/gsl/2.3/include
override LDFLAGS += -L/n/app/openblas/0.2.19/lib -L/n/app/gsl/2.3/lib/
endif

ND = nicksrc
NLIB = $(ND)/libnick.a

PROGS= qp3Pop qpDstat qpF4ratio qpAdm qpWave rolloff convertf dowtjack expfit.sh grabpars qpBound qpGraph qpreroot qpff3base rolloffp qpDpart

PROGS2 = rexpfit.r  wtjack.pl 

all: $(NLIB) $(PROGS) 

$(NLIB):
	$(MAKE) -C $(ND)

statsubs.o:     tables nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -o statsubs.o nicksrc/statsubs.c

expfit.sh:	expfit.templ
	sed -e "s&BBBIN&$(BIN)&" < expfit.templ > expfit.sh
	chmod +x expfit.sh

qpDstat: $(NLIB) qpDstat.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpDpart: $(NLIB) qpDpart.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qp3Pop: $(NLIB) qp3Pop.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpF4ratio: $(NLIB) qpF4ratio.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

rolloff: $(NLIB) rolloff.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

rolloffp: $(NLIB) rolloffp.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpBound: $(NLIB) qpBound.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o

dowtjack: $(NLIB) dowtjack.o

convertf: $(NLIB) convertf.o mcio.o admutils.o egsubs.o h2d.o exclude.o

grabpars: $(NLIB) grabpars.o

qpff3base: $(NLIB) qpff3base.o qpgsubs.o qpsubs.o mcio.o admutils.o egsubs.o regsubs.o ldsubs.o

qpAdm: $(NLIB) qpAdm.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o f4rank.o  eigsubs.o eigx.o

qpWave: $(NLIB) qpWave.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o f4rank.o eigsubs.o eigx.o

install:	all dirs tables
	mkdir -p $(TOP)
	cp $(PROGS) $(TOP)
	cp $(PROGS2) $(TOP)

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
	rm -f $(PROGS)

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/admixtables 

## gsl-config may help here


qpGraph: $(NLIB) qpGraph.o gslqp.o qpgsubs.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpreroot: qpreroot.o qpgsubs.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o
