CC = gcc

CFLAGS = -g -Wall 

LDFLAGS = -lm

EXECS = poisseq

all: $(EXECS)

poisseq: mainseq.o jacobiseq.o poisson.h
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o $(EXECS) $(TESTS) TAGS tags
