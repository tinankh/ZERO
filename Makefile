EXEC=./zero
SRCDIR=src
CFLAGS=-std=gnu99 -O3
Lib=-lpng -ltiff -ljpeg -lm $(OpenMP)
SRC_FILES=$(wildcard $(SRCDIR)/*.c)

all: $(EXEC)

$(EXEC): $(SRC_FILES)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

openmp: OpenMP=-fopenmp
openmp: zero

libzero: $(SRCDIR)/zero.c
	$(CC) -shared $^ -o $@.so -lm -fPIC -O3 -fopenmp

test: zero
	@echo
	@echo test on roma.png
	@echo ----------------
	./zero roma.png
	@echo
	@echo test on pelican.png
	@echo -------------------
	./zero pelican.png
	@echo
	@echo test on tampered1.png
	@echo ---------------------
	./zero tampered1.png
	@echo
	@echo test on tampered2.png
	@echo ---------------------
	./zero tampered2.png

clean:
	rm -f zero
	rm -f votes.png forgery.png forgery_c.png

.PHONY: clean
