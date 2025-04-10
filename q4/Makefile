CC = mpicc

CFLAGS = -g #-Wall -W #-D_CC_OVERLAP
LDFLAGS = -lm
POISSOBJS = decomp1d.o jacobi.o

EXECS = poiss2d  poiss2d_sr poiss2d_isr

all: $(EXECS)

poiss2d_isr: main_isr.o $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

poiss2d_sr: main_sr.o $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

poiss2d: main.o $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o $(EXECS) $(TESTS) TAGS tags
