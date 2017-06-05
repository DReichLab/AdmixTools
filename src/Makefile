HOMEL=$(PWD)
TOP=../bin
BIN=$(HOMEL)/../bin

override LDLIBS += -lgsl -lopenblas -lgfortran -lm -lnick
override LDFLAGS += -L./nicksrc
override CFLAGS += -c -g -p -Wimplicit -I./ -I./nicksrc
#CFLAGS= -c -g -p -Wimplicit-int

ND = nicksrc
NLIB = $(ND)/libnick.a

PROGS= qp3Pop qpDstat qpF4ratio qpAdm qpWave rolloff convertf dowtjack expfit.sh grabpars qpBound qpGraph
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

qpdslow: $(NLIB) qpdslow.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpdweight: $(NLIB) qpdweight.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpwest: $(NLIB) qpwest.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qp3Pop: $(NLIB) qp3Pop.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpF4ratio: $(NLIB) qpF4ratio.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

rolloff: $(NLIB) rolloff.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

qpBound: $(NLIB) qpBound.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o

dowtjack: $(NLIB) dowtjack.o

convertf: $(NLIB) convertf.o mcio.o admutils.o egsubs.o h2d.o exclude.o

grabpars: $(NLIB) grabpars.o
	
contam1: contam1.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o

mksamp1: $(NLIB) mksamp1.o qpsubs.o mcio.o admutils.o  egsubs.o h2d.o exclude.o regsubs.o

qpff3base: $(NLIB) qpff3base.o qpgsubs.o qpsubs.o mcio.o admutils.o egsubs.o regsubs.o ldsubs.o

qpginv: $(NLIB) qpginv.o qpgsubs.o qpsubs.o mcio.o admutils.o egsubs.o regsubs.o ldsubs.o 

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
