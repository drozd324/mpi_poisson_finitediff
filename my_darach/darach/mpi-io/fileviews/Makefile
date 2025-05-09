CC = gcc

MPICC=mpicc

CFLAGS = -g -Wall 

LDFLAGS =

EXECS = fvw1 fvr1 fvieww1 fviewr1

all: $(EXECS)

fvw1: fv1write.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

fvr1: fv1read.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

fvieww1: fvieww1.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

fviewr1: fviewr1.o
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)


.PHONY: clean

clean:
	$(RM) *.o $(EXECS)


