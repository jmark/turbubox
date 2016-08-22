.PHONY: modules

ifeq ($(findstring cheops,$(shell hostname)),cheops)

modules:
	@echo "Loading modules"
	$(call module load hdf5/1.8.13)
	$(call module load openmpi)
	$(call module load python/3.4.3 2> /dev/null)

else

modules: ;

endif

shell_avg.c: shell_avg.pyx modules
	cython $<

shell_avg.o: shell_avg.c modules
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
		-I/opt/rrzk/software/python/Python-3.4.3/include/python3.4m \
		-o $@ $<
