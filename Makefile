EXEC=./zero
SRCDIR=src
CFLAGS=-std=gnu99 -O3 $(OpenMP)
LIB=-lpng -ltiff -ljpeg -lm $(OpenMP)
SRC_FILES=$(wildcard $(SRCDIR)/*.c)

all: $(EXEC)

$(EXEC): $(SRC_FILES)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

openmp: OpenMP=-fopenmp
openmp: zero

test: zero
	@echo
	@echo test on roma.png
	@echo ----------------
	./zero roma.png
	@mv forgery.png roma_forgery.png
	@echo
	@echo test on pelican.png
	@echo -------------------
	./zero pelican.png
	@mv forgery.png pelican_forgery.png
	@echo
	@echo test on tampered1.png
	@echo ---------------------
	./zero tampered1.png
	@mv forgery.png tampered1_forgery.png
	@echo
	@echo test on tampered2.png
	@echo ---------------------
	./zero tampered2.png
	@mv forgery.png tampered2_forgery.png

clean:

	rm -f zero
	rm -f votes.png forgery.png forgery_c.png
	rm -f roma_forgery.png pelican_forgery.png
	rm -f tampered1_forgery.png tampered2_forgery.png

.PHONY: clean
