P4EST_INC = -I$(HOME)/opt/p4est/include
P4EST_LIB = -L$(HOME)/opt/p4est/lib -lp4est -lsc 
# P4WRAP= -L$(HOME)/projects/p4wrap/build -lp4wrap

CC		= gcc
CCFLAGS	= -O3 -std=c99 -pedantic -I.

# CCFLAGS += -DP4EST 
CCFLAGS += $(P4EST_INC)
CCFLAGS += $(P4EST_LIB)

all: libinterpolate.so libshellavg.so

libinterpolate.so: libinterpolate.c
	$(CC) $(CCFLAGS) -shared -fPIC -o $@ $< 

libshellavg.so: libshellavg.c
	$(CC) $(CCFLAGS) -shared -fPIC -o $@ $<

clean:
	rm -f libinterpolate.so
	rm -f libshellavg.so
