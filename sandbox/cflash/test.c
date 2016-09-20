#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "uniform.h"

int
main (int argc, char** argv)
{
	struct cflash_metadata_t meta;
	char* fpath;

	if (argc < 2) {
		fprintf(stderr, "missing file path\n");
		EXIT_CODE = 1;
		return EXIT_CODE;
	} else {
		fpath = argv[1];
	}

	cflash_fill_metadata(fpath,&meta);

	printf("%s -> %e\n", "gamma", meta.gamma);
	printf("%s -> %d\n", "nend", meta.nend);
	

	return EXIT_CODE;
}
