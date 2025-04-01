CC=mpicc
CFLAGS= -g -Wall
LDFLAGS=

EXECS = cartcreate1dshift

all: $(EXECS)

cartcreate1dshift: cartcreate1dshift.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) *.o $(EXECS)


