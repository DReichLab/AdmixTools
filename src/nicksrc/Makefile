override CFLAGS += -c -g -p -Wimplicit -fno-common -I./

all: libnick.a

libnick.a: strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o gds.o
	ar -r libnick.a $^
	ranlib libnick.a

clean:
	rm -f *.o
	rm -f libnick.a
