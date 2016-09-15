#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <string.h>

#define BUFLEN 80

int EXIT_CODE = 0;

union value_t {
	int i;
	float f;
	double d;
	char* s;
};

struct parameter_t {
	char name[BUFLEN+1];
	double value;
};

struct metadata_t {
	double gamma;
	double c_ambient;
	double rho_ambient;
	double domainsize[3];
	double min[3];
	double max[3];

	int nblock[3];

	int nend;
};

int
begins(char *str, char *key)
{
	return 0 == strncmp(str, key, strlen(key)) ? 1 : 0;
}

/* shuts whatever it gets thrown at and doesn't complain */
void
shut(int fd) {


}

void
transfer_parameter(char *name, union value_t value, struct metadata_t *meta)
{
		 if (begins(name,"gamma")) meta->gamma = value.d;
	else if (begins(name,"c_ambient")) meta->c_ambient = value.d;
	else if (begins(name,"rho_ambient")) meta->rho_ambient = value.d;

	else if (begins(name,"nend")) meta->nend = value.i;

	else if (begins(name,"bx")) meta->domainsize[0] = value.d;
	else if (begins(name,"by")) meta->domainsize[1] = value.d;
	else if (begins(name,"bz")) meta->domainsize[2] = value.d;

	else if (begins(name,"minx")) meta->min[0] = value.d;
	else if (begins(name,"miny")) meta->min[1] = value.d;
	else if (begins(name,"minz")) meta->min[2] = value.d;

	else if (begins(name,"maxx")) meta->max[0] = value.d;
	else if (begins(name,"maxy")) meta->max[1] = value.d;
	else if (begins(name,"maxz")) meta->max[2] = value.d;

	else if (begins(name,"nblockx")) meta->nblock[0] = value.i;
	else if (begins(name,"nblocky")) meta->nblock[1] = value.i;
	else if (begins(name,"nblockz")) meta->nblock[2] = value.i;
}

int
fill_parameter_list(hid_t dfile, char* dbname, struct metadata_t *meta)
{
	int retval = 0;
	hid_t dset = H5Dopen(dfile, dbname, H5P_DEFAULT);

	if (dset < 0) {retval=-1; goto clean_dset;}

	/* allocate array for data transfer */
	hsize_t dims[1];
	hid_t space = H5Dget_space(dset);
	H5Sget_simple_extent_dims(space, dims, NULL);
	struct parameter_t *rdata = (struct parameter_t*) malloc(dims[0] * sizeof (struct parameter_t));
	if (!rdata) {retval = -1; goto clean_rdata;}

	/* register name datatype */
	hid_t name_tid = H5Tcopy(H5T_C_S1);
	H5Tset_size(name_tid, BUFLEN*sizeof(char));

	/* get value datatype */
	hid_t cmpnd_tid = H5Dget_type(dset);
	hid_t value_tid = H5Tget_member_type(cmpnd_tid, 1);
	H5Tclose(cmpnd_tid);

	/* register compount data type as used in the FLASH file format */
	hid_t compound_tid = H5Tcreate(H5T_COMPOUND, sizeof(struct parameter_t));
	H5Tinsert(compound_tid,  "name", HOFFSET(struct parameter_t,  name),  name_tid);
	H5Tinsert(compound_tid, "value", HOFFSET(struct parameter_t, value), value_tid);

	/* transfer data into memory */
	H5Dread(dset, compound_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);

	H5Tclose(name_tid);
	H5Tclose(value_tid);
	H5Tclose(compound_tid);

	/* transfer data from bulky compount data array to slick struct */
	for (int i = 0; i <= dims[0]; i++)
		transfer_parameter(rdata[i].name, (union value_t) rdata[i].value, meta);

	clean_rdata:
		free(rdata);

	clean_dset:
		H5Dclose(dset);

	return retval;
}	


int
main (int argc, char** argv)
{
	hid_t dfile;
	struct metadata_t meta;
	char* fpath;

	if (argc < 2) {
		fprintf(stderr, "missing file path\n");
		EXIT_CODE = 1;
		return EXIT_CODE;
	} else {
		fpath = argv[1];
	}

	dfile = H5Fopen(fpath, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (dfile < 0) {EXIT_CODE=-1; goto clean_dfile;}

	fill_parameter_list(dfile, "integer runtime parameters", &meta);		
	fill_parameter_list(dfile, "real runtime parameters", &meta);		

	printf("%s -> %e\n", "gamma", meta.gamma);
	printf("%s -> %d\n", "nend", meta.nend);
	
	clean_dfile:
		H5Fclose(dfile);

	return EXIT_CODE;
}
