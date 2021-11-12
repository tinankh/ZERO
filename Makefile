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
	./zero roma.png roma99.jpg
	@echo
	@echo test on pelican.png
	@echo -------------------
	./zero pelican.png pelican99.jpg
	@echo
	@echo test on tampered1.png
	@echo ---------------------
	./zero tampered1.png tampered1_99.jpg
	@echo
	@echo test on tampered2.png
	@echo ---------------------
	./zero tampered2.png tampered2_99.jpg

clean:
	rm -f zero
	rm -f luminance.png votes.png votes_jpeg.png mask_f.png mask_m.png

.PHONY: clean
