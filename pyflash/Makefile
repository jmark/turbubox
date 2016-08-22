.PHONY: modules

ifeq ($(findstring cheops,$(shell hostname)),cheops)

modules:
	module load hdf5/1.8.13
	module load openmpi
	module load python/3.4.3 2> /dev/null

else

modules: ;

endif

shell_avg.c: modules shell_avg.pyx
	cython $<

shell_avg.so: modules shell_avg.c
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
		-I/opt/rrzk/software/python/Python-3.4.3/include/python3.4m \
		-o $@ $<

