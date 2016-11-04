CC		= gcc
CFLAGS	= -O3

libinterpolate.so: libinterpolate.c
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $<

libshellavg.so: libshellavg.c
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $<
