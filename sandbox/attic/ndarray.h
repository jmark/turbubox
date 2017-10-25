// dead simple ndarray implementation

typedef struct {
    int nshape;
    int *shape;
    int size;
    double *chunk;
    int ownchunk;
} nd_array;

int
nd_offset(const nd_array *const arr, const int *const indices, const int level)
{
    const int dim   = arr->shape[level];
    int index = indices[level];

    if (level < 0)
        return 0;

    // periodic indexing
    if (index < 0)
        index += dim;
    if (index >= dim)
        index -= dim;

    // note: arr[((i * ny) + j) * nz + k]
    return nd_offset(arr, indices, level-1) * dim + index;
}

//int
//nd_offset(const nd_array *const arr, const int *const indices)
//{
//    int acc = 0;
//    for (int i = arr->nshape-1 ; i >= 0 ; i--) {
//        const int *dim  = arr->shape[i];
//        const int index = indices[level];
//
//        if (level < 0)
//            return 0;
//
//        // periodic indexing
//        if (index < 0)
//            index += dim;
//        if (index >= dim)
//            index -= dim;
//
//        // note: arr[((i * ny) + j) * nz + k]
//        acc += * dim + index;
//    }
//}

inline double*
nd_getptr(const nd_array *const arr, const int *const indices)
{
    return arr->chunk + nd_offset(arr, indices, arr->nshape-1);
}

inline double
nd_get(const nd_array *const arr, const int *const indices)
{
    return *(nd_getptr(arr, indices));
}

inline void
nd_set(nd_array *const arr, const int *const indices, const double value)
{
    double *pt = nd_getptr(arr, indices);
    *pt = value;
}

inline void
nd_wipe(nd_array *arr) {
    if (arr->ownchunk)
        free(arr->chunk);
    free(arr->shape);
    free(arr);
}

inline nd_array*
nd_alloc(const int nshape, const int *const shape)
{
    nd_array *const arr = malloc(sizeof(nd_array));
    arr->shape    = malloc(sizeof(int) * nshape);
    arr->nshape   = nshape;

    int size = 1;
    for (int i = 0; i < nshape; i++) {
        arr->shape[i] = shape[i];
        size *= arr->shape[i];
    }

    arr->size = size;
    arr->chunk = 0;
    arr->ownchunk = 0;
    return arr;
}

inline nd_array*
nd_alloc_const(const int nshape, const int *const shape, const double konst)
{
    nd_array *arr = nd_alloc(nshape, shape);

    double *const chunk = malloc(sizeof(double) * arr->size);

    for (int i = 0; i < arr->size; i++)
        chunk[i] = konst;

    arr->chunk = chunk;
    arr->ownchunk = 1;

    return arr;
}

// ========================================================================= //

inline double*
nd_getptr_2(const nd_array *const arr, const int i, const int j)
{
    const int indices[] = {i,j};
    return nd_getptr(arr, indices);
}

inline double
nd_get_2(const nd_array *const arr, const int i, const int j)
{
    const int indices[] = {i,j};
    return nd_get(arr, indices);
}

inline void
nd_set_2(nd_array *const arr, const int i, const int j, const double value)
{
    const int indices[] = {i,j};
    nd_set(arr, indices, value);
}

inline nd_array*
nd_alloc_2(const int Nx, const int Ny)
{
    const int shape[] = {Nx, Ny};
    return nd_alloc(2, shape); 
}

inline nd_array*
nd_alloc_const_2(const int Nx, const int Ny, const double konst)
{
    const int shape[] = {Nx, Ny};
    return nd_alloc_const(2, shape, konst);
}

// ========================================================================= //

inline double*
nd_getptr_3(const nd_array *const arr, const int i, const int j, const int k)
{
    const int indices[] = {i,j,k};
    return nd_getptr(arr, indices);
}

inline double
nd_get_3(const nd_array *const arr, const int i, const int j, const int k)
{
    const int indices[] = {i,j,k};
    return nd_get(arr, indices);
}

inline void
nd_set_3(nd_array *const arr, const int i, const int j, const int k, const double value)
{
    const int indices[] = {i,j,k};
    nd_set(arr, indices, value);
}

inline nd_array*
nd_alloc_3(const int Nx, const int Ny, const int Nz)
{
    const int shape[] = {Nx, Ny, Nz};
    return nd_alloc(3, shape); 
}

// ========================================================================= //

inline double*
nd_getptr_4(const nd_array *const arr, const int i, const int j, const int k, const int l)
{
    const int indices[] = {i,j,k,l};
    return nd_getptr(arr, indices);
}

inline double
nd_get_4(const nd_array *const arr, const int i, const int j, const int k, const int l)
{
    const int indices[] = {i,j,k,l};
    return nd_get(arr, indices);
}

inline void
nd_set_4(nd_array *const arr, const int i, const int j, const int k, const int l, const double value)
{
    const int indices[] = {i,j,k,l};
    nd_set(arr, indices, value);
}

inline nd_array*
nd_alloc_4(const int Nx, const int Ny, const int Nz, const int Na)
{
    const int shape[] = {Nx, Ny, Nz, Na};
    return nd_alloc(4, shape); 
}
