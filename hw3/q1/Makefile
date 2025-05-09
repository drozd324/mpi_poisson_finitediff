CC = mpicc

CFLAGS =
DEBUG = -g -Wextra -Wall# -fsanitize=address -lefence
LDFLAGS = -lm
POISSOBJS = tools.o jacobi.o

EXECS = poiss2d #poiss2d_sr poiss2d_isr

all: $(EXECS)

#poiss2d_isr: main_isr.o $(POISSOBJS) 
#	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)
#
#poiss2d_sr: main_sr.o $(POISSOBJS) 
#	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)  
#
poiss2d: main.o $(POISSOBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)

.PHONY: clean tags tests

clean:
	rm -rf *.o $(EXECS)
