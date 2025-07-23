CC = gcc
CFLAGS = -Wall -pedantic -std=c99 -fopenmp -fcommon -g -O2
LDFLAGS = -lpthread -lm -ldl -lcblas -llapacke -lgsl -g -O2

EXDIR = bin
ODIR = src/obj

_OBJ = main.o eqm.o calibrate.o solver.o gnewton.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

model: $(OBJ)
	$(CC) -o $(EXDIR)/$@.exe $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f $(EXDIR)/*
