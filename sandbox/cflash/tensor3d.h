#define PAYLOAD_T double

struct t3d_obj_t {
	int dim[3];
	int len;

	PAYLOAD_T *vector;
};

t3d_obj_t*
t3d_alloc(int Nx,int Ny, int Nz)
{
	int len = Nx*Ny*Nz;

	t3d_obj_t *obj = (t3d_obj_t*) malloc(sizeof(t3d_obj_t));
	PAYLOAD_T *vec = (PAYLOAD_T*) malloc(len*sizeof(PAYLOAD_T));

	obj->len = len;

	obj->dim[0] = Nx;
	obj->dim[1] = Ny;
	obj->dim[2] = Nz;

	obj->vector = vec;

	return obj;
}

PAYLOAD_T*
t3d_getptr(t3d_obj_t* obj, int i, int j, int k)
{	
	int Nx = obj->dim[0];
	int Ny = obj->dim[1];
	int Nz = obj->dim[2];
	PAYLOAD_T *vec = obj->vector;

	//if (i >= Nx || j >= Ny || k >= k) raise error;

	return vec + (k*Ny + j)*Nx + i * sizeof(PAYLOAD_T);
}

PAYLOAD_T
t3d_get(t3d_obj_t* obj, int i, int j, int k)
{
	PAYLOAD_T *out = t3d_getptr(obj,i,j,k);
	return *out;
}

void
t3d_set(t3d_obj_t* obj, int i, int j, int k, PAYLOAD_T x)
{
	PAYLOAD_T *ptr = t3d_getptr(obj,i,j,k);
	*ptr = x;
}

void
t3d_plug_in(t3d_obj_t* self, t3d_obj_t* other, int I, int J, int K)
{
	int Nx = other->dim[0];
	int Ny = other->dim[1];
	int Nz = other->dim[2];


	PAYLOAD_T *vec = obj->vector;

	for (int i = I; i < Nx; i++)
	for (int j = J; j < Ny; j++)
	for (int k = K; k < Nx; k++) {
	}

	PAYLOAD_T *ptr = t3d_getptr(obj,i,j,k);
	*ptr = x;
}

void
t3d_scale(t3d_obj_t* obj, PAYLOAD_T scalar)
{
	int Nx = obj->dim[0];
	int Ny = obj->dim[1];
	int Nz = obj->dim[2];
	PAYLOAD_T *vec = obj->vector;

	for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nx; k++) {
		*vec *= scalar;
	}
}
