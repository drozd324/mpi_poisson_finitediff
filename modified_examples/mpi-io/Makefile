CC = gcc

MPICC=mpicc

CFLAGS = -g -Wall 

LDFLAGS =

EXECS = geni gend indpi1 indpd1 indpd2 indpdw1 readbini

all: $(EXECS)

geni: geni.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

gend: gend.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

indpi1: indpointersi.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

indpd1: indpointersd.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

indpd2: indpointersd2.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

indpdw1: indpointerswd1.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

readbini: readbini.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)


.PHONY: clean

clean:
	$(RM) *.o $(EXECS)


