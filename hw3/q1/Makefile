MPICC = mpicc
CC = $(MPICC)

CFLAGS =
DEBUG = -Wextra -Wall # -fsanitize=address -lefence
LDFLAGS = -lm

EXECS = poiss2d 

all: $(EXECS)

poiss2d: main.c tools.o jacobi.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.o $(EXECS)
