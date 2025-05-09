CC=mpicc

CFLAGS= -g -Wall

LDFLAGS=

EXECS = cartcreate1dshift  cartcreate1dget cartcreate1dcoords #cartcreate1dshiftlr

all: $(EXECS)


cartcreate1dshift: cartcreate1dshift.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

cartcreate1dshiftlr: cartcreate1dshiftlr.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

cartcreate1dget: cartcreate1dget.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

cartcreate1dcoords: cartcreate1dcoords.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) *.o $(EXECS)


