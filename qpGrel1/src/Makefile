HOMEL=$(PWD)
TOP=../bin
## binary for install  
# DEBUG_OPTIONS= -g -pg
DEBUG_OPTIONS= -g -pg
BIN=$(HOMEL)/../bin
PERLDIR=$(HOMEL)/perlsrc

NLIB=$(HOMEL)/admixlib/nicklib.a
IDIR=$(HOMEL)/admixinclude
VPATH=.:nicksrc
DEBUG_OPTIONS= -g -pg

CFLAGS= -c -g -p -I$(IDIR) -Wimplicit

OBJ=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o
TWTAB=\"$(HOMEL)/admixtables/twtable\"

statsubs.o:     tables	nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o nicksrc/statsubs.c


T2=dowtjack
T20=dowtjack.o

M1=qpDstat 
M1O=qpDstat.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o
 
M2=qp3Pop 
M2O=qp3Pop.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o
 
M3=qpF4ratio  
M3O=qpF4ratio.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o
 
M4=rolloff  
M4O=rolloff.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o
 
M5=qpBound  
M5O=qpBound.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o
 
X1=convertf  
X1O=convertf.o  mcio.o  admutils.o  egsubs.o h2d.o  exclude.o

X2=grabpars  
X2O=grabpars.o 

expfit.sh:	expfit.templ
	sed -e "s&BBBIN&$(BIN)&" < expfit.templ > expfit.sh
	chmod +x expfit.sh
 
IDIR=$(HOMEL)/admixinclude  
##ISTR= -I$(HOME)/include  -I/opt/nag/libC/include
##NAGL= -L/opt/nag/libC/lib -lnagc_nag -lm -pthread

##nagqp.o:		nagqp.c
#	$(CC)  $(CFLAGS) -I/opt/nag/libC/include -o nagqp.o nagqp.c


$(M1): nicklib	$(M1O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M1) $(M1O) $(NLIB)


$(T2): nicklib $(T20)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(T2) $(T20) $(NLIB)



$(M2): nicklib	$(M2O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M2) $(M2O) $(NLIB)

$(M3): nicklib	$(M3O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M3) $(M3O) $(NLIB)

$(M4): nicklib	$(M4O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M4) $(M4O) $(NLIB)

$(M5): nicklib	$(M5O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(M5) $(M5O) $(NLIB)

$(X1): nicklib	$(X1O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(X1) $(X1O) $(NLIB)

$(X2): nicklib	$(X2O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(X2) $(X2O) $(NLIB)

#NAGL= -L/opt/nag/libC/lib -lnagc_nag -lm -pthread

PROGS= qp3Pop qpDstat qpF4ratio rolloff convertf dowtjack expfit.sh grabpars qpBound
PROGS2 = rexpfit.r  wtjack.pl 
PROGS3 = qpGraph qpreroot

all: nicklib all3 $(PROGS) 
all3: $(PROGS3) 

install:	all	tables
	mkdir -p $(TOP)
	cp $(PROGS) $(TOP)
	cp $(PROGS2) $(TOP)

uninstall:  
	rm -f $(NLIB) tables
	rm -f $(TOP)/*

libnick.a:	dirs tables  $(OBJ)
	ar -r libnick.a $(OBJ)
	ranlib libnick.a

nicklib:	dirs  tables  libnick.a 
	cp libnick.a  $(NLIB)

tables:    
	echo "tables made"  > tables
	cp twtable  $(HOMEL)/admixtables
	
dirs:	
	mkdir -p  $(HOMEL)/admixlib
	mkdir -p  $(HOMEL)/admixtables
	mkdir -p  $(HOMEL)/admixinclude
	mkdir -p  $(BIN)
	cp  *.h  $(IDIR)
	cp  nicksrc/*.h  $(IDIR)

clean: 
	rm -f *.o 
	rm -f *junk*
	rm -f core
	rm -f libnick.a
	rm -f $(PROGS)
	rm -f nicksrc/*.o

clobber: clean rmdirs uninstall


rmdirs: 
	rm -rf $(HOMEL)/admixlib 
	rm -rf $(HOMEL)/admixtables 
	rm -rf $(HOMEL)/admixinclude 

ISTR= -I$(HOME)/include  -I/opt/nag/libC/include
NAGL= -L/opt/nag/libC/lib -lnagc_nag -lm -pthread
CFLAGS= -c -g -p $(ISTR) -Wimplicit-int

M9F=qpGraph 
M9FO=qpGraph.o nagqp.o qpgsubs.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o

M9G=qpreroot 
M9GO=qpreroot.o qpgsubs.o qpsubs.o mcio.o ldsubs.o admutils.o  egsubs.o regsubs.o

$(M9F): nicklib $(M9FO)
	gcc -I$(IDIR) $(NAGL) $(DEBUG_OPTIONS) -lm  -o $(M9F) $(M9FO) $(NLIB)

$(M9G): nicklib $(M9GO)
	gcc -I$(IDIR)  $(DEBUG_OPTIONS) -lm  -o $(M9G) $(M9GO) $(NLIB)

