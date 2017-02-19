# include <stdio.h>
# include <math.h>
# include <malloc.h>

double
LagrangePolynomial(const double *xs, const int xslen, const int j, const double X) {
    double acc = 1;
    
    for (int i = 0; i < xslen; ++i) {
        if (i==j) continue;
        acc *= (X-xs[i])/(xs[j]-xs[i]);
    }

    return acc;
}

void
lagrange_interpolate_2d_RG(
    const int xslen, const double *xs, const double *fs,
    const int Xslen, const double *Xs,       double *Fs
)
{
    double *Ls = malloc(sizeof(double) * Xslen * xslen);

    for (int I = 0; I < Xslen; I++)
    for (int i = 0; i < xslen; i++)
        Ls[I*xslen + i] = LagrangePolynomial(xs, xslen, i, Xs[I]);

    for (int I = 0; I < Xslen; I++)
    for (int J = 0; J < Xslen; J++) {
        double F = 0;

        for (int i = 0; i < xslen; i++)
        for (int j = 0; j < xslen; j++)
            F += fs[(i * xslen) + j] * Ls[I*xslen + i]*Ls[J*xslen + j];

        Fs[(I * Xslen) + J] = F;
    }
    free(Ls);
}

void
lagrange_interpolate_3d_RG(
    const int xslen, const double *xs, const double *fs,
    const int Xslen, const double *Xs,       double *Fs
)
{
    double *Ls = malloc(sizeof(double) * Xslen * xslen);

    for (int I = 0; I < Xslen; I++)
    for (int i = 0; i < xslen; i++)
        Ls[I*xslen + i] = LagrangePolynomial(xs, xslen, i, Xs[I]);

    for (int I = 0; I < Xslen; I++)
    for (int J = 0; J < Xslen; J++)
    for (int K = 0; K < Xslen; K++) {
        double F = 0;

        for (int i = 0; i < xslen; i++)
        for (int j = 0; j < xslen; j++)
        for (int k = 0; k < xslen; k++)
            F += fs[((i * xslen) + j) * xslen + k] * Ls[I*xslen + i]*Ls[J*xslen + j]*Ls[K*xslen + k];

        Fs[((I * Xslen) + J) * Xslen + K] = F;
    }
    free(Ls);
}

void
box_to_flexi(
    const int xslen, const double *xs, 
    const int Nx, const int Ny, const int Nz, const double *fss,
    const int Xslen, const double *Xs, double *Fss,
    const int oflen, const int *offsets
)
{
    // init Lagrange Interpolation Matrix
    double *Ls = malloc(sizeof(double) * Xslen * xslen);

    for (int I = 0; I < Xslen; I++)
    for (int i = 0; i < xslen; i++)
        Ls[I*xslen + i] = LagrangePolynomial(xs, xslen, i, Xs[I]);

    const int Offset = Xslen * Xslen * Xslen;

    // loop over elements provided by FLEXI
    for (int elemid = 0; elemid < oflen; elemid++) {
        const double *fs = fss + offsets[elemid]; // non-consecutive
              double *Fs = Fss + Offset *elemid;  // consecutive

        // interpolate
        for (int I = 0; I < Xslen; I++)
        for (int J = 0; J < Xslen; J++)
        for (int K = 0; K < Xslen; K++) {
            double F = 0;

            for (int i = 0; i < xslen; i++)
            for (int j = 0; j < xslen; j++)
            for (int k = 0; k < xslen; k++) {
                //const double f = fs[((i * xslen) + j) * xslen + k];
                const double f = fs[((i * Ny) + j) * Nz + k];
                F += f * Ls[I*xslen + i]*Ls[J*xslen + j]*Ls[K*xslen + k];
                //printf("%d %d %d -> %f\n", i,j,k, f);
            }
            Fs[((I * Xslen) + J) * Xslen + K] = F;
        }
    }

    free(Ls);
}

# define PIX(index, dim) (\
    (index) < 0      ? (index) + (dim) : \
    (index) >= (dim) ? (index) - (dim) : index)

void
box_to_elements(
    const int Nx, const int Ny, int Nz, double *box, 
    int nelems, double *indices, int nx, int ny, int nz, double *elems, int offset)
{
    const int stride = nx*ny*nz;
    for (int elemid = 0; elemid < nelems; elemid++) {

        const double *const index = indices + elemid * 3;
        const int I = round(index[0]);
        const int J = round(index[1]);
        const int K = round(index[2]);

        double *const elem = elems + elemid * stride;

        for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
            const int _i = PIX(I + i - offset, Nx);
            const int _j = PIX(J + j - offset, Ny);
            const int _k = PIX(K + k - offset, Nz);

            elem[((i * ny) + j) * nz + k] = box[((_i * Ny) + _j) * Nz + _k];
        }
    }
}

void
elements_to_box(
    const int Nx, const int Ny, int Nz, double *box,
    int nelems, double *indices, int nx, int ny, int nz, double *elems)
{
    const int stride = nx*ny*nz;
    for (int elemid = 0; elemid < nelems; elemid++) {

        const double *const index = indices + elemid * 3;
        const int I = round(index[0]);
        const int J = round(index[1]);
        const int K = round(index[2]);

        double *const elem = elems + elemid * stride;

        for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
            const int _i = I + i;
            const int _j = J + j;
            const int _k = K + k;

            box[((_i * Ny) + _j) * Nz + _k] = elem[((i * ny) + j) * nz + k];
        }
    }
}

void
elements_to_box_fv(
    const int Nx, const int Ny, int Nz, double *box,
    int nelems, double *indices, int nx, int ny, int nz, double *elems,
    int *fvs)
{
    const int stride = nx*ny*nz;
    for (int elemid = 0; elemid < nelems; elemid++) {

        if (fvs[elemid] == 0)
            continue;

        const double *const index = indices + elemid * 3;
        const int I = round(index[0]);
        const int J = round(index[1]);
        const int K = round(index[2]);

        double *const elem = elems + elemid * stride;

        for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
            const int _i = I + i;
            const int _j = J + j;
            const int _k = K + k;

            box[((_i * Ny) + _j) * Nz + _k] = elem[((i * ny) + j) * nz + k];
        }
    }
}

void
box_to_elements_avg_boundaries(
    const int Nx, const int Ny, int Nz, double *box, 
    int nelems, double *indices, int nx, int ny, int nz, double *elems)
{
    const int stride = nx*ny*nz;
    for (int elemid = 0; elemid < nelems; elemid++) {

        const double *const index = indices + elemid * 3;
        const int I = round(index[0]);
        const int J = round(index[1]);
        const int K = round(index[2]);

        double *const elem = elems + elemid * stride;

        for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {

            const int _i = I + i;
            const int _j = J + j;
            const int _k = K + k;

            double acc = box[((_i * Ny) + _j) * Nz + _k];
            double cnt = 1.0;

            if (i == 0) {
                acc += box[((PIX(_i-1,Nx) * Ny) + _j) * Nz + _k];
                cnt++;
            }
            if (j == 0) {
                acc += box[((_i * Ny) + PIX(_j-1,Ny)) * Nz + _k];
                cnt++;
            }
            if (k == 0) {
                acc += box[((_i * Ny) + _j) * Nz + PIX(_k-1,Nz)];
                cnt++;
            }

            if (i == nx-1) {
                acc += box[((PIX(_i+1,Nx) * Ny) + _j) * Nz + _k];
                cnt++;
            }
            if (j == ny-1) {
                acc += box[((_i * Ny) + PIX(_j+1,Ny)) * Nz + _k];
                cnt++;
            }
            if (k == nz-1) {
                acc += box[((_i * Ny) + _j) * Nz + PIX(_k+1,Nz)];
                cnt++;
            }

            elem[((i * ny) + j) * nz + k] = acc/cnt;
        }
    }
}

void
change_grid_space(
    const int nelems,
    const int nx, const int ny, const int nz ,const double *xs, double *fss, 
    const int Nx, const int Ny, const int Nz ,const double *Xs, double *Fss)
{
    double *Ls = malloc(sizeof(double) * Nx * nx);
    for (int I = 0; I < Nx; I++)
    for (int i = 0; i < nx; i++)
        Ls[I*nx + i] = LagrangePolynomial(xs, nx, i, Xs[I]);

    const int stride = nx*ny*nz;
    const int Stride = Nx*Ny*Nz;

    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *const fs = fss + elemid * stride;
              double *const Fs = Fss + elemid * Stride;

        for (int I = 0; I < Nx; I++)
        for (int J = 0; J < Ny; J++)
        for (int K = 0; K < Nz; K++) {
            double F = 0;

            for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                F += fs[((i * ny) + j) * nz + k] * Ls[I*nx + i]*Ls[J*ny + j]*Ls[K*nz + k];

            Fs[((I * Ny) + J) * Nz + K] = F;
        }
    }
}

void
change_grid_space_fv(
    const int nelems,
    const int nx, const int ny, const int nz ,const double *xs, double *fss, 
    const int Nx, const int Ny, const int Nz ,const double *Xs, double *Fss,
    const int *fvs)
{
    double *Ls = malloc(sizeof(double) * Nx * nx);
    for (int I = 0; I < Nx; I++)
    for (int i = 0; i < nx; i++)
        Ls[I*nx + i] = LagrangePolynomial(xs, nx, i, Xs[I]);

    const int stride = nx*ny*nz;
    const int Stride = Nx*Ny*Nz;

    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *const fs = fss + elemid * stride;
              double *const Fs = Fss + elemid * Stride;

        //printf("elemid,fv: %d, %d\n", elemid, fvs[elemid]);
        for (int I = 0; I < Nx; I++)
        for (int J = 0; J < Ny; J++)
        for (int K = 0; K < Nz; K++) {
            double F = 0;

            if (fvs[elemid] > 0) {
                F = fs[((I * Ny) + J) * Nz + K];
            } else {
                for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    F += fs[((i * ny) + j) * nz + k] * Ls[I*nx + i]*Ls[J*ny + j]*Ls[K*nz + k];
            }

            Fs[((I * Ny) + J) * Nz + K] = F;
        }
    }
}



void
flexi_to_box(
    const int xslen, const double *xs,
    const int Xslen, const double *Xs, 
    const int nelems, const int *offsets, const double *flexi,
    const int Nx, const int Ny, const int Nz, double *box
)
{
    // init Lagrange Interpolation Matrix
    double *Ls = malloc(sizeof(double) * Xslen * xslen);

    for (int I = 0; I < Xslen; I++)
    for (int i = 0; i < xslen; i++)
        Ls[I*xslen + i] = LagrangePolynomial(xs, xslen, i, Xs[I]);

    const int stride = xslen * xslen * xslen;

    // loop over elements provided by FLEXI
    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *fs = flexi + stride *elemid;
              double *Fs = box   + offsets[elemid];

        // interpolate
        for (int I = 0; I < Xslen; I++)
        for (int J = 0; J < Xslen; J++)
        for (int K = 0; K < Xslen; K++) {
            double F = 0;

            for (int i = 0; i < xslen; i++)
            for (int j = 0; j < xslen; j++)
            for (int k = 0; k < xslen; k++) {
                //printf("foo: ");
                const double f = fs[((i * xslen) + j) * xslen + k];
                F += f * Ls[I*xslen + i]*Ls[J*xslen + j]*Ls[K*xslen + k];
                //printf("%d %d %d -> %f\n", i,j,k, f);
            }
            //printf("%d / %d -> %d\n", elemid, nelems, ((I * Ny) + J) * Nz + K);
            Fs[((I * Ny) + J) * Nz + K] = F;
            //printf("%d %d %d -> %f\n", I,J,K, F);
        }
    }

    free(Ls);
}
