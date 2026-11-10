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

# 4 regions (USA, CHN, low-tariff, high-tariff partners): bin/model_nc4.bin, reads
# ../python/output/iomat_nc4.txt (icio.py --nc4); outputs get _nc4
ODIRNC4 = src/obj_nc4
OBJNC4 = $(patsubst %,$(ODIRNC4)/%,$(_OBJ))

$(ODIRNC4)/%.o: src/%.c $(DEPS)
	@mkdir -p $(ODIRNC4)
	$(CC) -c -o $@ $< $(CFLAGS) -DNC=4

model_nc4: $(OBJNC4)
	$(CC) -o $(EXDIR)/$@.bin $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o $(ODIR200)/*.o $(ODIRNC4)/*.o
	rm -f $(EXDIR)/*
