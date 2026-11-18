CC = gcc
CFLAGS = -Wall -pedantic -std=c99 -fopenmp -fcommon -g
LDFLAGS = -lpthread -lm -ldl -lcblas -llapacke -lgsl -g

EXDIR = bin
ODIR = src/obj
DEPS = $(wildcard src/*.h)

_OBJ = main.o eqm.o calibrate.o solver.o gnewton.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# bin/model.bin: the baseline model (100-period horizon); bin/model_nt200.bin: 200-period
# horizon, used for the portfolio-adjustment-cost runs (-b), outputs get _nt200
all: model model_nt200

$(ODIR)/%.o: src/%.c $(DEPS)
	@mkdir -p $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS)

model: $(OBJ)
	@mkdir -p $(EXDIR) output
	$(CC) -o $(EXDIR)/$@.bin $^ $(CFLAGS) $(LDFLAGS)

ODIR200 = src/obj200
OBJ200 = $(patsubst %,$(ODIR200)/%,$(_OBJ))

$(ODIR200)/%.o: src/%.c $(DEPS)
	@mkdir -p $(ODIR200)
	$(CC) -c -o $@ $< $(CFLAGS) -DNT=200

model_nt200: $(OBJ200)
	@mkdir -p $(EXDIR) output
	$(CC) -o $(EXDIR)/$@.bin $^ $(CFLAGS) $(LDFLAGS)

.PHONY: all clean

clean:
	rm -f $(ODIR)/*.o $(ODIR200)/*.o
	rm -f $(EXDIR)/*
