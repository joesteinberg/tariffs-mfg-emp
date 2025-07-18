CC = gcc
CFLAGS = -O2 -Wall -pedantic -std=c99 -fopenmp -fcommon -g
LDFLAGS = -lpthread -lm -ldl -lcblas -llapacke -lgsl -g

EXDIR = bin
ODIR = src/obj

_OBJ = main.o eqm.o calibrate.o exporters.o solver.o gnewton.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

model: $(OBJ)
	$(CC) -o $(EXDIR)/$@.exe $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f $(EXDIR)/*
