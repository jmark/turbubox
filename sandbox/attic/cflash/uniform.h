/*
 * _cflash_* ... private variables/functions
 * cflash_*  ... public api
 *
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <string.h>

#define _CFLASH_BUFLEN 80

int EXIT_CODE = 0;

union _cflash_value_t {
	int i;
	float f;
	double d;
	char* s;
};

struct cflash_metadata_t {
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
_cflash_begins(char *str, char *key)
{
	return 0 == strncmp(str, key, strlen(key)) ? 1 : 0;
}

void
_cflash_transfer_parameter(char *name, union _cflash_value_t value, struct cflash_metadata_t *meta)
{
		 if (_cflash_begins(name,"gamma")) meta->gamma = value.d;
	else if (_cflash_begins(name,"c_ambient")) meta->c_ambient = value.d;
	else if (_cflash_begins(name,"rho_ambient")) meta->rho_ambient = value.d;

	else if (_cflash_begins(name,"nend")) meta->nend = value.i;

	else if (_cflash_begins(name,"bx")) meta->domainsize[0] = value.d;
	else if (_cflash_begins(name,"by")) meta->domainsize[1] = value.d;
	else if (_cflash_begins(name,"bz")) meta->domainsize[2] = value.d;

	else if (_cflash_begins(name,"minx")) meta->min[0] = value.d;
	else if (_cflash_begins(name,"miny")) meta->min[1] = value.d;
	else if (_cflash_begins(name,"minz")) meta->min[2] = value.d;

	else if (_cflash_begins(name,"maxx")) meta->max[0] = value.d;
	else if (_cflash_begins(name,"maxy")) meta->max[1] = value.d;
	else if (_cflash_begins(name,"maxz")) meta->max[2] = value.d;

	else if (_cflash_begins(name,"nblockx")) meta->nblock[0] = value.i;
	else if (_cflash_begins(name,"nblocky")) meta->nblock[1] = value.i;
	else if (_cflash_begins(name,"nblockz")) meta->nblock[2] = value.i;
}

int
_cflash_fill_parameter_list(hid_t dfile, char* dbname, struct cflash_metadata_t *meta)
{
	struct parameter_t {
		char name[_CFLASH_BUFLEN+1];
		double value;
	};

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
	H5Tset_size(name_tid, _CFLASH_BUFLEN*sizeof(char));

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
		_cflash_transfer_parameter(rdata[i].name, (union _cflash_value_t) rdata[i].value, meta);

	clean_rdata:
		free(rdata);

	clean_dset:
		H5Dclose(dset);

	return retval;
}	

void
cflash_fill_metadata(char *fpath, struct cflash_metadata_t *meta)
{
	hid_t dfile;

	dfile = H5Fopen(fpath, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (dfile < 0) {EXIT_CODE=-1; goto clean_dfile;}

	_cflash_fill_parameter_list(dfile, "integer runtime parameters", meta);		
	_cflash_fill_parameter_list(dfile, "real runtime parameters", meta);		

	clean_dfile:
		H5Fclose(dfile);
}
