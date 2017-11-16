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
change_basis_3d(
    const int nelems, const int nn,
    const double *Vdm, const double *fss, double *Fss)
{
    const int stride = nn*nn*nn;

    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *const fs = fss + elemid * stride;
              double *const Fs = Fss + elemid * stride;

        for (int I = 0; I < nn; I++)
        for (int J = 0; J < nn; J++)
        for (int K = 0; K < nn; K++) {
            double F = 0.0;

            for (int i = 0; i < nn; i++)
            for (int j = 0; j < nn; j++)
            for (int k = 0; k < nn; k++)
                F += fs[(i*nn + j) * nn + k] * Vdm[I*nn + i]*Vdm[J*nn + j]*Vdm[K*nn + k];
                
            Fs[((I * nn) + J) * nn + K] = F;
        }
    }
}

void
change_basis_3d_2(
    const int nelems, const int nn,
    const double *Vdm, const double *fss, double *Fss)
{
    const int stride = nn*nn*nn;

    double  *xi = malloc(sizeof(double) * stride);
    double *eta = malloc(sizeof(double) * stride);
    double *psi = malloc(sizeof(double) * stride);

    // loop over elements
    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *const fs = fss + elemid * stride;
              double *const Fs = Fss + elemid * stride;

        // nullify temporary arrays
        for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
        for (int k = 0; k < nn; k++)
            xi[(i*nn + j) * nn + k] = eta[(i*nn + j) * nn + k] = psi[(i*nn + j) * nn + k] = 0.0;

        // xi direction
        for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
        for (int k = 0; k < nn; k++)
        for (int l = 0; l < nn; l++)
             xi[(i*nn + j) * nn + k] +=  fs[(l*nn + j) * nn + k] * Vdm[i*nn + l];

        // eta direction
        for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
        for (int k = 0; k < nn; k++)
        for (int l = 0; l < nn; l++)
            eta[(i*nn + j) * nn + k] +=  xi[(i*nn + l) * nn + k] * Vdm[j*nn + l];

        // psi direction
        for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
        for (int k = 0; k < nn; k++)
        for (int l = 0; l < nn; l++)
            psi[(i*nn + j) * nn + k] += eta[(i*nn + j) * nn + l] * Vdm[k*nn + l];

        // copy to result array
        for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
        for (int k = 0; k < nn; k++)
            Fs[((i * nn) + j) * nn + k] = psi[(i*nn + j) * nn + k];
    }

    free(xi); free(eta); free(psi);
}

void
change_grid_space_2d(
    const int nelems,
    const int nx, const int ny, const double *xs, const double *fss, 
    const int Nx, const int Ny, const double *Xs, double *Fss)
{
    double *Ls = malloc(sizeof(double) * Nx * nx);
    for (int I = 0; I < Nx; I++)
    for (int i = 0; i < nx; i++)
        Ls[I*nx + i] = LagrangePolynomial(xs, nx, i, Xs[I]);

    const int stride = nx*ny;
    const int Stride = Nx*Ny;

    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *const fs = fss + elemid * stride;
              double *const Fs = Fss + elemid * Stride;

        for (int I = 0; I < Nx; I++)
        for (int J = 0; J < Ny; J++) {
            double F = 0;

            for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                F += fs[i*ny + j] * Ls[I*nx + i]*Ls[J*ny + j];

            Fs[(I * Ny) + J] = F;
        }
    }
    free(Ls);
}

void
change_grid_space(
    const int nelems,
    const int nx, const int ny, const int nz ,const double *xs, const double *fss, 
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
                F += fs[(i*ny + j) * nz + k] * Ls[I*nx + i]*Ls[J*ny + j]*Ls[K*nz + k];

            Fs[((I * Ny) + J) * Nz + K] = F;
        }
    }
    free(Ls);
}

void
change_grid_space_dg_fv(
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

            if (fvs[elemid] > 0) { // fv element
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
    free(Ls);
}

void
change_grid_space_fv_dg(
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

            if (fvs[elemid] == 0) { // dg element
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
    free(Ls);
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

# define I1(nx,i)                     (i)
# define I2(nx,ny,i,j)               ((i)*(ny)) + (j)
# define I3(nx,ny,nz,i,j,k)         (((i)*(ny)  + (j))*(nz) + (k))
# define I4(nx,ny,nz,nw,i,j,k,l)   ((((i)*(ny)  + (j))*(nz) + (k)))*(nw) + (l)

void
blocks_to_box(
    const int  rlevel,  const int nblocks,
    const int *rlevels, const double *coords, const double *domain,
    const int nx, const int ny, const int nz, const double *blocks,
    const int Nx, const int Ny, const int Nz, double *box
) {
    double domsize[6];
    double pos[3]; // position vector
    int    iul[3]; // index vector upper left 

    domsize[I1(3,0)] = domain[I2(2,3,1,0)] - domain[I2(2,3,0,0)];
    domsize[I1(3,1)] = domain[I2(2,3,1,1)] - domain[I2(2,3,0,1)];
    domsize[I1(3,2)] = domain[I2(2,3,1,2)] - domain[I2(2,3,0,2)];

    for (int rl = 0; rl < nblocks; rl++) {

        if (rlevel != rlevels[rl]) continue;

        pos[0] = ((double)Nx)*(coords[I2(nblocks,3,rl,0)] - domain[I2(2,3,0,0)])/domsize[I1(3,0)];
        pos[1] = ((double)Ny)*(coords[I2(nblocks,3,rl,1)] - domain[I2(2,3,0,1)])/domsize[I1(3,1)];
        pos[2] = ((double)Nz)*(coords[I2(nblocks,3,rl,2)] - domain[I2(2,3,0,2)])/domsize[I1(3,2)];

        iul[0] = round(pos[0]) - nx/2;
        iul[1] = round(pos[1]) - ny/2;
        iul[2] = round(pos[2]) - nz/2;

        // fill box with matching blocks
        for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            box[I3(Nx,Ny,Nz,iul[0]+i,iul[1]+j,iul[2]+k)] = blocks[I4(nblocks,nx,ny,nz,rl,i,j,k)];
    }
}
