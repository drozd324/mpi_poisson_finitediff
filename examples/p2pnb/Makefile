CC = mpicc

#MPE_ROOT=/home/shared/darach/install/mpe2
MPE_ROOT=

#PROFLAGS= -I$(MPE_ROOT)/include -L$(MPE_ROOT)/lib

CFLAGS = -g -Wall 

LIBDIR=$(MPE_ROOT)/lib

#PROFLIBS=$(LIBDIR)/mpe_prof.o $(LIBDIR)/libmpe.a -D_PROFMPE

LDFLAGS = # -lmpe

EXECS =  nbsr1

all: $(EXECS)

nbsr1: basicnbsendrecv.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) *.o $(EXECS)


