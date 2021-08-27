CFLAGS= -O3
LIB=-lpng -ltiff -ljpeg -lm $(OpenMP)
SRC_FILES=$(wildcard src/*.c)

all: zero

zero: $(SRC_FILES)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

openmp: OpenMP=-fopenmp
openmp: zero

libzero: src/zero.c src/iio.c
	$(CC) -shared $^ -o $@.so -lm -fPIC -O3 -fopenmp $(LIB)

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
	rm -f luminance.png votes.png forgery.png forgery_c.png main_grid.txt

.PHONY: clean
