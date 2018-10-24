include ../compilers.inc

CC=cc
CFLAGS=-O3 -hfp3 -I${HDF5_DIR}/include
LDFLAGS=-L${HDF5_DIR}/lib
LIBS=-lhdf5

cleanall:
	make clean
	make all

all: h5_convert file_compare

h5_convert: h5_convert.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o h5_convert h5_convert.c $(LIBS)

file_compare: file_compare.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o file_compare file_compare.c $(LIBS)

clean:
	rm file_compare h5_convert
