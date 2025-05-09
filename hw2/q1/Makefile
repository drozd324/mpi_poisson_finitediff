CC = mpicc

CFLAGS = -g -Wall -W #-D_CC_OVERLAP
LDFLAGS = -lm
POISSOBJS = decomp1d.o jacobi.o

EXECS = poiss1d #poiss1dnb

all: $(EXECS)

poiss1d: main.o $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

poiss1dnb: mainnb.o $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o $(EXECS) $(TESTS) TAGS tags
