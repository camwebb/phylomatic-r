all: phylomatic.so

phylomatic.so:
	R CMD SHLIB -O -o phylomatic.so phylomatic.c

debug:
	MAKEFLAGS="CFLAGS=-g\ -O1\ -march=x86-64\ -mtune=generic\ -pipe\ -fno-plt" R CMD SHLIB -o phylomatic.so phylomatic.c

clean:
	rm -f *.o *.so
