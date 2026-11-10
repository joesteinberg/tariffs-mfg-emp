CC = gcc
CFLAGS = -Wall -pedantic -std=c99 -fopenmp -fcommon -g
LDFLAGS = -lpthread -lm -ldl -lcblas -llapacke -lgsl -g

EXDIR = bin
ODIR = src/obj
DEPS = $(wildcard src/*.h)

_OBJ = main.o eqm.o calibrate.o solver.o gnewton.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

model: $(OBJ)
	$(CC) -o $(EXDIR)/$@.bin $^ $(CFLAGS) $(LDFLAGS)

# 200-period horizon (portfolio-adjustment-cost runs, -b): bin/model_nt200.bin, outputs get _nt200
ODIR200 = src/obj200
OBJ200 = $(patsubst %,$(ODIR200)/%,$(_OBJ))

$(ODIR200)/%.o: src/%.c $(DEPS)
	@mkdir -p $(ODIR200)
	$(CC) -c -o $@ $< $(CFLAGS) -DNT=200

model_nt200: $(OBJ200)
	$(CC) -o $(EXDIR)/$@.bin $^ $(CFLAGS) $(LDFLAGS)

# 3 regions (USA, CHN, rest of world; the original version): bin/model_nc3.bin, reads
# ../python/output/iomat_nc3.txt (icio.py); outputs get _nc3
ODIRNC3 = src/obj_nc3
OBJNC3 = $(patsubst %,$(ODIRNC3)/%,$(_OBJ))

$(ODIRNC3)/%.o: src/%.c $(DEPS)
	@mkdir -p $(ODIRNC3)
	$(CC) -c -o $@ $< $(CFLAGS) -DNC=3

model_nc3: $(OBJNC3)
	$(CC) -o $(EXDIR)/$@.bin $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(ODIR200)/*.o $(ODIRNC3)/*.o
	rm -f $(EXDIR)/*
