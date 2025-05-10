MPICC=mpicc
CFLAGS =  -Wall -Wextra 
LDFLAGS = 
EXECS = main 

all: $(EXECS)

main: main.c
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean

clean:
	rm -rf *.o $(EXECS)


