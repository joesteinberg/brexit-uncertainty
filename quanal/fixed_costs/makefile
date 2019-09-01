CC = gcc
CFLAGS = -g -Wall -pedantic -std=c99 -fopenmp -m64 -I/opt/intel/mkl/lib/intel64 -I/opt/intel/mkl/include
LDFLAGS = -lgsl -lcblas -L/opt/intel/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl

EXDIR = bin
ODIR = src/obj

_OBJ = main.o eqm.o calibrate.o exporters.o solver.o gnewton.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

brexit: $(OBJ)
	$(CC) -o $(EXDIR)/$@ $^ $(CFLAGS) $(LDFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f $(EXDIR)/*
