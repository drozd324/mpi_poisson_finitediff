CC=mpicc
CFLAGS= -Wall -W
LDFLAGS=

EXECS = cartcreate1dshift  cartcreate2dshift

all: $(EXECS)

cartcreate1dshift: cartcreate1dshift.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

cartcreate2dshift: cartcreate2dshift.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean

clean:
	rm -r *.o $(EXECS)


