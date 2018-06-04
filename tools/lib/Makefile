P4EST       = -DP4EST -I$(HOME)/opt/include -L$(HOME)/opt/lib -lp4est -lsc 
P4WRAP      = -L$(HOME)/projects/p4wrap/build -lp4wrap

CC		= gcc
CFLAGS	= -O3 -std=c99 -pedantic -I.

all: libinterpolate.so libshellavg.so

libinterpolate.so: libinterpolate.c
	$(CC) $(CFLAGS) $(P4EST) -shared -fPIC -o $@ $<

libshellavg.so: libshellavg.c
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $<

clean:
	rm -f libinterpolate.so
	rm -f libshellavg.so
