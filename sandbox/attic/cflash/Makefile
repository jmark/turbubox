CC = gcc
CCFLAGS = -std=c99 -Wall
LDFLAGS = -lhdf5

test: test.c uniform.h
	$(CC) $(CCFLAGS) $(LDFLAGS) -o $@ $<
