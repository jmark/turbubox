#P4EST       = -DP4EST -I$(HOME)/opt/p4est/include -L$(HOME)/opt/p4est/lib -lp4est -lsc 
#P4WRAP      = -L$(HOME)/projects/p4wrap/build -lp4wrap

CC		= gcc
CCFLAGS	= -O3 -std=c99 -pedantic -I.

# CCFLAGS += -DP4EST 
# CCFLAGS += -I$(P4EST_INC)
# CCFLAGS += -L$(P4EST_LIB)
# CCFLAGS += -lp4est -lsc

all: libinterpolate.so libshellavg.so

libinterpolate.so: libinterpolate.c
	$(CC) $(CCFLAGS) -shared -fPIC -o $@ $< 

libshellavg.so: libshellavg.c
	$(CC) $(CCFLAGS) -shared -fPIC -o $@ $<

clean:
	rm -f libinterpolate.so
	rm -f libshellavg.so
