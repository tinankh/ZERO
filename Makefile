IIO= -std=c99 -lpng -ltiff -ljpeg
OPT= -O3
OpenMP=
FSAN=

zero: zero.c iio.c iio.h
	$(CC) $(OPT) $(OpenMP) -o $@ zero.c iio.c $(IIO) -lm $(FSAN)

openmp: OpenMP=-fopenmp
openmp: zero

sanitize: OpenMP=-fopenmp
sanitize: FSAN=-fsanitize=thread -g
sanitize: zero

test: zero
	@echo
	@echo test on roma.pgm
	@echo ----------------
	./zero roma.pgm
	@mv forgery.png roma_forgery.png
	@echo
	@echo test on pelican.ppm
	@echo -------------------
	./zero pelican.ppm
	@mv forgery.png pelican_forgery.png
	@echo
	@echo test on tampered1.pgm
	@echo ---------------------
	./zero tampered1.pgm
	@mv forgery.png tampered1_forgery.png
	@echo
	@echo test on tampered2.ppm
	@echo ---------------------
	./zero tampered2.ppm
	@mv forgery.png tampered2_forgery.png

clean:
	rm -f zero
	rm -f votes.png forgery.png forgery_c.png
	rm -f roma_forgery.png pelican_forgery.png
	rm -f tampered1_forgery.png tampered2_forgery.png
