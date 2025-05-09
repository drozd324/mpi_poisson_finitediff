CC = mpicc

CFLAGS = -g -Wall

LDFLAGS = -lm

EXECS = vectortype1 vectortype2 general1only general1only2 general1only_bcast \
	generalarray generalarraywrong mpigetcountvsgetelements general1p

all: $(EXECS)

general1only: general1only.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)
general1only2: general1only2.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)
general1only_bcast: general1only_bcast.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

general1p: general1p.c
	$(CC) $(CFLAGS) -std=c99 -o $@ $< $(LDFLAGS)

generalarray: generalarray.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)
generalarraywrong: generalarraywrong.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

vectortype1: vectortype1.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

vectortype2: vectortype2.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

contig1: contig1.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

mpigetcountvsgetelements: mpigetcountvsgetelements.o
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean

clean:
	$(RM) *.o $(EXECS)


