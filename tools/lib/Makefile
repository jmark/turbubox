CC		= gcc
CCFLAGS	= -O3 -std=c99 -pedantic -I.

CCFLAGS += -I$(P4EST_INC)
CCFLAGS += -L$(P4EST_LIB)

all: libinterpolate.so libshellavg.so

libinterpolate.so: libinterpolate.c
	$(CC) -DP4EST $(CCFLAGS) -shared -fPIC -o $@ $< -lp4est -lsc

libshellavg.so: libshellavg.c
	$(CC) $(CCFLAGS) -shared -fPIC -o $@ $<

clean:
	rm -f libinterpolate.so
	rm -f libshellavg.so
