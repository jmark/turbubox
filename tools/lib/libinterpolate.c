# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include <stdint.h>

# if defined(P4EST)
# include <p4est.h>
# include <p4est_connectivity.h>

# include <p8est.h>
# include <p8est_connectivity.h>
# endif

# define I1(nx,i)                     (i)
# define I2(nx,ny,i,j)               ((i)*(ny)) + (j)
# define I3(nx,ny,nz,i,j,k)         (((i)*(ny)  + (j))*(nz) + (k))
# define I4(nx,ny,nz,nw,i,j,k,l)   ((((i)*(ny)  + (j))*(nz) + (k)))*(nw) + (l)

# define PI 3.1415926535897932384626433832795

/* ========================================================================= */

double inline dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double inline norm3(const double a[3], const double b[3])
{
    return sqrt(dot3(a,b));
}

double inline dist3(const double a[3], const double b[3])
{
    const double c0 = b[0]-a[0];
    const double c1 = b[1]-a[1];
    const double c2 = b[2]-a[2];

    return sqrt(c0*c0 + c1*c1 + c2*c2);
}

void cross3(const double a[3], const double b[3], double *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

void pclose3(const double p[3], const double n[3], const double x[3], double *y)
{
    /* 
     * p[3] ... positition vector of plane
     * n[2] ... normal vector of plane
     * x[3] ... point nearby plane
     * y[3] ... closest point on plane
     */

    //printf("dot3: %f\n", dot3(n,n));
    const double fac = (dot3(x,n) - dot3(p,n))/dot3(n,n);

    y[0] = x[0] - fac*n[0];
    y[1] = x[1] - fac*n[1];
    y[2] = x[2] - fac*n[2];
}

int inline imin(const int a, const int b)
{
    return a < b ? a : b;
}

int inline imax(const int a, const int b)
{
    return a > b ? a : b;
}

/* ========================================================================= */

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
change_grid_space_2d_2(
    const int nelems,
    const int nx, const int ny,
    const int Nx, const int Ny,
    const double *Lss, const double *fss, double *Fss)
{
    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *const fs = fss + elemid * nx*ny;
              double *const Fs = Fss + elemid * Nx*Ny;

        for (int I = 0; I < Nx; I++)
        for (int J = 0; J < Ny; J++) {
            double F = 0;

            for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                F += fs[i*ny + j] * Lss[I4(Nx,Ny,nx,ny,I,J,i,j)];
                
            Fs[(I * Ny) + J] = F;
        }
    }
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

    for (int elemid = 0; elemid < nelems; elemid++) {
        const double *const fs = fss + elemid * nx*ny;
              double *const Fs = Fss + elemid * Nx*Ny;

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

inline int linearsearch(const int len, const double nodes[], const double x)
{
    for (int i = 1; i < len; i++)
        if (nodes[i] > x)
            return i-1;

    return len-1;
}

void narrow_down(
    const int nx, const double xnodes[], const double x, int *iL, int *iR)
{
    const int L = linearsearch(nx, xnodes, x);
    const int R = L+1 >= nx ? L : L+1;

    *iL = L;
    *iR = R;
}

inline double nearest2D(
    const int nx, const double xnodes[],
    const int ny, const double ynodes[],
    const double fs[], const double x, const double y)
{
    int ix = linearsearch(nx, xnodes, x);
    int iy = linearsearch(ny, ynodes, y);

    if (ix+1 < nx)
        ix = fabs(xnodes[ix]-x) < fabs(xnodes[ix+1]-x) ? ix : ix+1;

    if (iy+1 < ny)
        iy = fabs(ynodes[iy]-y) < fabs(ynodes[iy+1]-y) ? iy : iy+1;

    return fs[I2(nx,ny,ix,iy)];
}

inline double nearest3D(
    const int nx, const double xnodes[],
    const int ny, const double ynodes[],
    const int nz, const double znodes[],
    const double fs[], const double x, const double y, const double z)
{
    int ix = linearsearch(nx, xnodes, x);
    int iy = linearsearch(ny, ynodes, y);
    int iz = linearsearch(nz, znodes, z);

    if (ix+1 < nx)
        ix = fabs(xnodes[ix]-x) < fabs(xnodes[ix+1]-x) ? ix : ix+1;

    if (iy+1 < ny)
        iy = fabs(ynodes[iy]-y) < fabs(ynodes[iy+1]-y) ? iy : iy+1;

    if (iz+1 < nz)
        iz = fabs(znodes[iz]-z) < fabs(znodes[iz+1]-z) ? iz : iz+1;

    return fs[I3(nx,ny,nz,ix,iy,iz)];
}

inline double linear_interpolation(
    const double x1, const double x2,
    const double f1, const double f2,
    const double x
){
    return (f1*(x2-x) + f2*(x-x1))/(x2-x1);
}

inline double bilinear_interpolation(
    const double x1,  const double x2,
    const double y1,  const double y2, 
    const double f11, const double f21,
    const double f12, const double f22,
    const double x,   const double y
){
    return (y2-y)/(y2-y1) * ((x2-x)/(x2-x1)*f11 + (x-x1)/(x2-x1)*f21)
         + (y-y1)/(y2-y1) * ((x2-x)/(x2-x1)*f12 + (x-x1)/(x2-x1)*f22);
}

inline double trilinear_interpolation(
    const double x1,  const double x2,
    const double y1,  const double y2, 
    const double z1,  const double z2, 
    const double f111, const double f211,
    const double f121, const double f221,
    const double f112, const double f212,
    const double f122, const double f222,
    const double x,   const double y,   const double z
){
    const double f1 = bilinear_interpolation(x1,x2,y1,y2,f111,f211,f121,f221,x,y);
    const double f2 = bilinear_interpolation(x1,x2,y1,y2,f112,f212,f122,f222,x,y);
    return linear_interpolation(z1,z2,f1,f2,z);
}

inline double bicosine_interpolation(
    const double x1,  const double x2,
    const double y1,  const double y2, 
    const double f11, const double f21,
    const double f12, const double f22,
    const double x,   const double y
){
    const double X = x1 + 0.5*(x2-x1)*(1 - cos(PI*(x-x1)/(x2-x1)));
    const double Y = y1 + 0.5*(y2-y1)*(1 - cos(PI*(y-y1)/(y2-y1)));

    return (y2-Y)/(y2-y1) * ((x2-X)/(x2-x1)*f11 + (X-x1)/(x2-x1)*f21)
         + (Y-y1)/(y2-y1) * ((x2-X)/(x2-x1)*f12 + (X-x1)/(x2-x1)*f22);
}

inline double bilinear(
    const int nx, const double xnodes[],
    const int ny, const double ynodes[],
    const double fs[], const double x, const double y)
{
    const int i = linearsearch(nx,xnodes,x);
    const int j = linearsearch(ny,ynodes,y);

    /* Assume nx >= 2. */
    const int iL = i+2 < nx ? i : nx-2;
    const int iR = iL+1;

    /* Assume ny >= 2. */
    const int jL = j+2 < ny ? j : ny-2;
    const int jR = jL+1;

    const double x1 = xnodes[iL];
    const double x2 = xnodes[iR];

    const double y1 = ynodes[jL];
    const double y2 = ynodes[jR];

    const double f11 = fs[I2(nx,ny,iL,jL)];
    const double f21 = fs[I2(nx,ny,iR,jL)];
    const double f12 = fs[I2(nx,ny,iL,jR)];
    const double f22 = fs[I2(nx,ny,iR,jR)];

    return bilinear_interpolation(x1,x2,y1,y2,f11,f21,f12,f22,x,y);
}

// inline double bilinear(
//     const int nx, const double xnodes[],
//     const int ny, const double ynodes[],
//     const double fs[], const double x, const double y)
// {
//     const int ix = linearsearch(nx, xnodes, x);
//     const int iy = linearsearch(ny, ynodes, y);
// 
//     const int ox = nx < 2 ? 0 : 1;
//     const int oy = ny < 2 ? 0 : 1;
// 
//     const double x1 = nx < 2 ? 0 : xnodes[ix];
//     const double x2 = nx < 2 ? 1 : xnodes[ix+1];
// 
//     const double y1 = ny < 2 ? 0 : ynodes[iy];
//     const double y2 = ny < 2 ? 1 : ynodes[iy+1];
// 
//     const double f11 = fs[I2(nx,ny,ix,   iy)];
//     const double f21 = fs[I2(nx,ny,ix+ox,iy)];
//     const double f12 = fs[I2(nx,ny,ix,   iy+oy)];
//     const double f22 = fs[I2(nx,ny,ix+ox,iy+oy)];
// 
//     return bilinear_interpolation(x1,x2,y1,y2,f11,f21,f12,f22,x,y);
// }

inline double trilinear(
    const int nx, const double xnodes[],
    const int ny, const double ynodes[],
    const int nz, const double znodes[],
    const double fs[], const double x, const double y, const double z)
{
    const int i = linearsearch(nx,xnodes,x);
    const int j = linearsearch(ny,ynodes,y);
    const int k = linearsearch(nz,znodes,z);

    /* Assume nx >= 2. */
    const int iL = i+2 < nx ? i : nx-2;
    const int iR = iL+1;

    /* Assume ny >= 2. */
    const int jL = j+2 < ny ? j : ny-2;
    const int jR = jL+1;

    /* Assume nz >= 2. */
    const int kL = k+2 < nz ? k : nz-2;
    const int kR = kL+1;

    const double x1 = xnodes[iL];
    const double x2 = xnodes[iR];

    const double y1 = ynodes[jL];
    const double y2 = ynodes[jR];

    const double z1 = znodes[kL];
    const double z2 = znodes[kR];

    const double f111 = fs[I3(nx,ny,nz,iL,jL,kL)];
    const double f211 = fs[I3(nx,ny,nz,iR,jL,kL)];
    const double f121 = fs[I3(nx,ny,nz,iL,jR,kL)];
    const double f221 = fs[I3(nx,ny,nz,iR,jR,kL)];

    const double f112 = fs[I3(nx,ny,nz,iL,jL,kR)];
    const double f212 = fs[I3(nx,ny,nz,iR,jL,kR)];
    const double f122 = fs[I3(nx,ny,nz,iL,jR,kR)];
    const double f222 = fs[I3(nx,ny,nz,iR,jR,kR)];

    return trilinear_interpolation(x1,x2,y1,y2,z1,z2,f111,f211,f121,f221,f112,f212,f122,f222,x,y,z);
}

inline double bicosine(
    const int nx, const double xnodes[],
    const int ny, const double ynodes[],
    const double fs[], const double x, const double y)
{
    const int ix = linearsearch(nx, xnodes, x);
    const int iy = linearsearch(ny, ynodes, y);

    const int ox = nx < 2 ? 0 : 1;
    const int oy = ny < 2 ? 0 : 1;

    const double x1 = nx < 2 ? 0 : xnodes[ix];
    const double x2 = nx < 2 ? 1 : xnodes[ix+1];

    const double y1 = ny < 2 ? 0 : ynodes[iy];
    const double y2 = ny < 2 ? 1 : ynodes[iy+1];

    const double f11 = fs[I2(nx,ny,ix,   iy)];
    const double f21 = fs[I2(nx,ny,ix+ox,iy)];
    const double f12 = fs[I2(nx,ny,ix,   iy+oy)];
    const double f22 = fs[I2(nx,ny,ix+ox,iy+oy)];

    return bicosine_interpolation(x1,x2,y1,y2,f11,f21,f12,f22,x,y);
}

# define GRIDLINES 0

enum {
    NEAREST, LINEAR, COSINE
};

void
cells_to_image_2d(
    const int nnodetype[1], const int *nodetype,
    const int ncoords[2], const double *coords,
    const int nsizes[2], const double *sizes,
    const int ncells[4], const double *cells,
    const int nimage[2], double *const image,
    const int method
) {
    const int nc = ncells[0];
    const int nx = ncells[1];
    const int ny = ncells[2];

    const int idx = nimage[0];
    const int idy = nimage[1];

    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);

    for (size_t icell = 0; icell < nc; icell++) {

        // Render only LEAF nodes.
        if (nodetype[icell] != 1) continue;

        const double xlength = sizes[I2(nc,3,icell,0)];
        const double ylength = sizes[I2(nc,3,icell,1)];

        const double xvert = coords[I2(nc,3,icell,0)] - 0.5*xlength;
        const double yvert = coords[I2(nc,3,icell,1)] - 0.5*ylength;

        // printf("%ld %f %f %f %f\n",icell,xlength,ylength,xvert,yvert);

        for (int i = 0; i < nx; i++)
            xnodes[i] = xvert + (i+0.5)/nx * xlength;

        for (int i = 0; i < ny; i++)
            ynodes[i] = yvert + (i+0.5)/ny * ylength;

        const int imgx = idx *  xvert            + 0.5*xlength/nx;
        const int Imgx = idx * (xvert + xlength) - 0.5*xlength/nx + 1;

        const int imgy = idy *  yvert            + 0.5*ylength/ny;
        const int Imgy = idy * (yvert + ylength) - 0.5*ylength/ny + 1;

        for (int i = imgx; i < Imgx; i++)
        for (int j = imgy; j < Imgy; j++)
        {
            const double x = (i+0.5)/idx;
            const double y = (j+0.5)/idy;
            
            switch(method) {
                case NEAREST:
                    image[I2(idx,idy,i,j)] = nearest2D(nx, xnodes, ny,ynodes, &cells[I3(nc,nx,ny,icell,0,0)],x,y);
                    break;

                case LINEAR:
                    image[I2(idx,idy,i,j)] = bilinear(nx, xnodes, ny,ynodes, &cells[I3(nc,nx,ny,icell,0,0)],x,y);
                    break;

                case COSINE:
                    image[I2(idx,idy,i,j)] = bicosine(nx, xnodes, ny,ynodes, &cells[I3(nc,nx,ny,icell,0,0)],x,y);
                    break;
            }
        }
    }

    free(xnodes);
    free(ynodes);
}

void
cells_to_image_3d(
    const int ncoords[2], const double *coords,
    const int nsizes[2], const double *sizes,
    const int ncells[4], const double *cells,
    const int nimage[3], double *const image,
    const int method
) {
    const int nc = ncells[0];
    const int nx = ncells[1];
    const int ny = ncells[2];
    const int nz = ncells[3];

    const int idx = nimage[0];
    const int idy = nimage[1];
    const int idz = nimage[2];

    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);
    double *const znodes = malloc(sizeof(double) * nz);

    for (size_t icell = 0; icell < nc; icell++) {
        const double xlength = sizes[I2(nc,3,icell,0)];
        const double ylength = sizes[I2(nc,3,icell,1)];
        const double zlength = sizes[I2(nc,3,icell,2)];

        const double xvert = coords[I2(nc,3,icell,0)] - 0.5*xlength;
        const double yvert = coords[I2(nc,3,icell,1)] - 0.5*ylength;
        const double zvert = coords[I2(nc,3,icell,2)] - 0.5*zlength;

        //printf("%ld | %f %f %f | %f %f %f\n",icell,xlength,ylength,zlength,xvert,yvert,zvert);

        for (int i = 0; i < nx; i++)
            xnodes[i] = xvert + (i+0.5)/nx * xlength;

        for (int i = 0; i < ny; i++)
            ynodes[i] = yvert + (i+0.5)/ny * ylength;

        for (int i = 0; i < nz; i++)
            znodes[i] = zvert + (i+0.5)/nz * zlength;

        const int imgx = idx *  xvert            + 0.5*xlength/nx     + GRIDLINES;
        const int Imgx = idx * (xvert + xlength) - 0.5*xlength/nx + 1 - GRIDLINES;

        const int imgy = idy *  yvert            + 0.5*ylength/ny     + GRIDLINES;
        const int Imgy = idy * (yvert + ylength) - 0.5*ylength/ny + 1 - GRIDLINES;

        const int imgz = idz *  zvert            + 0.5*zlength/nz     + GRIDLINES;
        const int Imgz = idz * (zvert + zlength) - 0.5*zlength/nz + 1 - GRIDLINES;

        for (int i = imgx; i < Imgx; i++)
        for (int j = imgy; j < Imgy; j++)
        for (int k = imgz; k < Imgz; k++)
        {
            const double x = (i+0.5)/idx;
            const double y = (j+0.5)/idy;
            const double z = (k+0.5)/idz;

            //image[I3(idx,idy,idz,i,j,k)] = nearest3D(
            //    nx,xnodes,ny,ynodes,nz,znodes,&cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);

            switch(method) {
                case NEAREST:
                    image[I3(idx,idy,idz,i,j,k)] = nearest3D(
                        nx,xnodes,ny,ynodes,nz,znodes,&cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);
                    break;

                case LINEAR:
                    image[I3(idx,idy,idz,i,j,k)] = trilinear(
                        nx,xnodes,ny,ynodes,nz,znodes,&cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);
                    break;
            }
        }
    }

    free(xnodes);
    free(ynodes);
    free(znodes);
}

void
morton_to_coords(
    const int ncenters[2], const double *centers,
    const int nsizes[2], const double *sizes,
    const int ncoords[2], double *coords
) {
    const int nc = ncoords[0]; // # of cells

    for (size_t icell = 0; icell < nc; icell++) {
        const double xlength = sizes[I2(nc,2,icell,0)];
        const double ylength = sizes[I2(nc,2,icell,1)];

        const double xvert = centers[I2(nc,2,icell,0)] - 0.5*xlength;
        const double yvert = centers[I2(nc,2,icell,1)] - 0.5*ylength;

        coords[I2(nc,3,icell,0)] = xvert;
        coords[I2(nc,3,icell,1)] = yvert;
        coords[I2(nc,3,icell,2)] = xlength;
    }
}

void
cells_to_plane_3d(
    const int nnodetype[1], const int *nodetype,
    const int ncoords[2], const double *coords,
    const int nsizes[2], const double *sizes,
    const int ncells[4], const double *cells,
    const int nimage[2], double *const image,
    const double p[3], const double u[3], const double v[3],
    const int method
) {
    const int nc = ncells[0];
    const int nx = ncells[1];
    const int ny = ncells[2];
    const int nz = ncells[3];

    // physical locations of nodes within cube
    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);
    double *const znodes = malloc(sizeof(double) * nz);

    const double uu = dot3(u,u);
    const double vv = dot3(v,v);

    const double su = sqrt(uu);
    const double sv = sqrt(vv);

    double normal[3];
    cross3(u,v,normal);

    for (size_t icell = 0; icell < nc; icell++) {

        // Render only LEAF nodes.
        if (nodetype[icell] != 1) continue;

        double pclose[3];
        double vshift[3];
        double center[3];

        const double xlength = sizes[I2(nc,3,icell,0)];
        const double ylength = sizes[I2(nc,3,icell,1)];
        const double zlength = sizes[I2(nc,3,icell,2)];

        const double xhalf = 0.5*xlength;
        const double yhalf = 0.5*ylength;
        const double zhalf = 0.5*zlength;

        const double xdelt = xlength/nx;
        const double ydelt = ylength/ny;
        const double zdelt = zlength/nz;

        const double xcenter = coords[I2(nc,3,icell,0)];
        const double ycenter = coords[I2(nc,3,icell,1)];
        const double zcenter = coords[I2(nc,3,icell,2)];

        const double xvert = xcenter - xhalf;
        const double yvert = ycenter - yhalf;
        const double zvert = zcenter - zhalf;

        /* ----------------------------------------------------------------- */

        center[0] = xcenter;
        center[1] = ycenter;
        center[2] = zcenter;

        pclose3(p,normal,center,pclose);

        if (dist3(center,pclose) > xlength)
            continue;
            
        /* ----------------------------------------------------------------- */

        for (int i = 0; i < nx; i++)
            xnodes[i] = xvert + (i+0.5)/nx * xlength;

        for (int i = 0; i < ny; i++)
            ynodes[i] = yvert + (i+0.5)/ny * ylength;

        for (int i = 0; i < nz; i++)
            znodes[i] = zvert + (i+0.5)/nz * zlength;

        /* ----------------------------------------------------------------- */

# if 0
        // vshift[0] = verts[0]-p[0];
        // vshift[1] = verts[1]-p[1];
        // vshift[2] = verts[2]-p[2];

        vshift[0] = pclose[0] - 2*length*(u[0]/su + v[0]/sv);
        vshift[1] = pclose[1] - 2*length*(u[1]/su + v[1]/sv);
        vshift[2] = pclose[2] - 2*length*(u[2]/su + v[2]/sv);

        const int ix = imax(0,dimage[0] * dot3(vshift,u)/uu - 0.5);
        const int iy = imax(0,dimage[1] * dot3(vshift,v)/vv - 0.5);

        /* ----------------------------------------------------------------- */

        // vshift[0] = verts[0]-p[0] + length;
        // vshift[1] = verts[1]-p[1] + length;
        // vshift[2] = verts[2]-p[2] + length;

        vshift[0] = pclose[0] + 2*length*(u[0]/su + v[0]/sv);
        vshift[1] = pclose[1] + 2*length*(u[1]/su + v[1]/sv);
        vshift[2] = pclose[2] + 2*length*(u[2]/su + v[2]/sv);

        const int Ix = imin(dimage[0]-1,dimage[0] * dot3(vshift,u)/uu - 0.5);
        const int Iy = imin(dimage[0]-1,dimage[1] * dot3(vshift,v)/vv - 0.5);

        //printf("%d %d %d %d\n", ix, Ix, iy, Iy);
# else
        const int ix = 0;
        const int iy = 0;

        const int Ix = nimage[0]-1;
        const int Iy = nimage[1]-1;
# endif

        /* ----------------------------------------------------------------- */

        for (int i = ix; i <= Ix; i++)
        for (int j = iy; j <= Iy; j++)
        {
            const double x = p[0] + (i+0.5)/nimage[0]*u[0] + (j+0.5)/nimage[1]*v[0];
            const double y = p[1] + (i+0.5)/nimage[0]*u[1] + (j+0.5)/nimage[1]*v[1];
            const double z = p[2] + (i+0.5)/nimage[0]*u[2] + (j+0.5)/nimage[1]*v[2];

            if (!(fabs(xcenter-x) <= xdelt+xhalf
               && fabs(ycenter-y) <= ydelt+yhalf
               && fabs(zcenter-z) <= zdelt+zhalf)) {
                continue;
            }

            switch(method) {
                case NEAREST:
                    image[I2(nimage[0],nimage[1],i,j)]
                        = nearest3D(nx,xnodes,ny,ynodes,nz,znodes,&cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);
                    break;

                case LINEAR:
                    image[I2(nimage[0],nimage[1],i,j)]
                        = trilinear(nx,xnodes,ny,ynodes,nz,znodes,&cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);
                    break;
            }
        }
    }

    free(xnodes);
    free(ynodes);
    free(znodes);
}

// cells_to_plane_3d(
//     const int ncoords[2], const double *coords,
//     const int nsizes[2], const double *sizes,
//     const int ncells[4], const double *cells,
//     const int nimage[2], double *const image,
//     const double p[3], const double u[3], const double v[3],
//     const int method
// )

int intersection(
    const double p1[3], const double p2[3], const double pE[3], const double nE[3], double *x
) {
    double u[3],q[3];

    u[0] = p2[0] - p1[0];
    u[1] = p2[1] - p1[1];
    u[2] = p2[2] - p1[2];

    const double t = dot3(u,nE);

    if (fabs(t) < 1e-9)
        return 0;

    q[0] = pE[0] - p1[0];
    q[1] = pE[1] - p1[1];
    q[2] = pE[2] - p1[2];

    const double r = dot3(q,nE)/t;

    if (-0.001 < r && r < 1.001) {
        x[0] = p1[0] + r * u[0];
        x[1] = p1[1] + r * u[1];
        x[2] = p1[2] + r * u[2];

        return 1;   
    }

    return 0;
}

int
plane_morton_to_coords(
    const int nnodetype[1], const int *nodetype,
    const int ncoords[2], const double *coords,
    const int nsizes[2], const double *sizes,
    const double p[3], const double u[3], const double v[3],
    const int nlines[3], double *lines, const int doedges
) {
    const int nc = ncoords[0]; // # of cells

    const double uu = dot3(u,u);
    const double vv = dot3(v,v);

    const double lu = sqrt(uu);
    const double lv = sqrt(vv);

    double n[3];
    cross3(u,v,n);

    int ecount = 0;

    for (size_t icell = 0; icell < nc; icell++) {

        // Render only LEAF nodes.
        if (nodetype[icell] != 1) continue;

        double verts[3];
        double center[3];
        double pclose[3];

        const double length = sizes[I2(nc,3,icell,0)];
        const double half = 0.5*length;

        double edges[6*4*2*3];

        /* ----------------------------------------------------------------- */

        center[0] = coords[I2(nc,3,icell,0)];
        center[1] = coords[I2(nc,3,icell,1)];
        center[2] = coords[I2(nc,3,icell,2)];

        verts[0] = center[0] - 0.5*sizes[I2(nc,3,icell,0)];
        verts[1] = center[1] - 0.5*sizes[I2(nc,3,icell,1)];
        verts[2] = center[2] - 0.5*sizes[I2(nc,3,icell,2)];

        pclose3(p,n,center,pclose);

        if (!(fabs(center[0]-pclose[0]) <= half
           && fabs(center[1]-pclose[1]) <= half
           && fabs(center[2]-pclose[2]) <= half)) {
            continue;
        }

        /* north */
        edges[I4(6,4,2,3, 0,0,0,0)] = verts[0];
        edges[I4(6,4,2,3, 0,0,0,1)] = verts[1];
        edges[I4(6,4,2,3, 0,0,0,2)] = verts[2];

        edges[I4(6,4,2,3, 0,0,1,0)] = verts[0];
        edges[I4(6,4,2,3, 0,0,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 0,0,1,2)] = verts[2];


        edges[I4(6,4,2,3, 0,1,0,0)] = verts[0];
        edges[I4(6,4,2,3, 0,1,0,1)] = verts[1];
        edges[I4(6,4,2,3, 0,1,0,2)] = verts[2];

        edges[I4(6,4,2,3, 0,1,1,0)] = verts[0];
        edges[I4(6,4,2,3, 0,1,1,1)] = verts[1];
        edges[I4(6,4,2,3, 0,1,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 0,2,0,0)] = verts[0];
        edges[I4(6,4,2,3, 0,2,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 0,2,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 0,2,1,0)] = verts[0];
        edges[I4(6,4,2,3, 0,2,1,1)] = verts[1];
        edges[I4(6,4,2,3, 0,2,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 0,3,0,0)] = verts[0];
        edges[I4(6,4,2,3, 0,3,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 0,3,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 0,3,1,0)] = verts[0];
        edges[I4(6,4,2,3, 0,3,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 0,3,1,2)] = verts[2];


        /* souTh */
        edges[I4(6,4,2,3, 1,0,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,0,0,1)] = verts[1];
        edges[I4(6,4,2,3, 1,0,0,2)] = verts[2];

        edges[I4(6,4,2,3, 1,0,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,0,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 1,0,1,2)] = verts[2];


        edges[I4(6,4,2,3, 1,1,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,1,0,1)] = verts[1];
        edges[I4(6,4,2,3, 1,1,0,2)] = verts[2];

        edges[I4(6,4,2,3, 1,1,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,1,1,1)] = verts[1];
        edges[I4(6,4,2,3, 1,1,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 1,2,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,2,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 1,2,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 1,2,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,2,1,1)] = verts[1];
        edges[I4(6,4,2,3, 1,2,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 1,3,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,3,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 1,3,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 1,3,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 1,3,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 1,3,1,2)] = verts[2];



        /* wesT */
        edges[I4(6,4,2,3, 2,0,0,0)] = verts[0];
        edges[I4(6,4,2,3, 2,0,0,1)] = verts[1];
        edges[I4(6,4,2,3, 2,0,0,2)] = verts[2];

        edges[I4(6,4,2,3, 2,0,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 2,0,1,1)] = verts[1];
        edges[I4(6,4,2,3, 2,0,1,2)] = verts[2];


        edges[I4(6,4,2,3, 2,1,0,0)] = verts[0];
        edges[I4(6,4,2,3, 2,1,0,1)] = verts[1];
        edges[I4(6,4,2,3, 2,1,0,2)] = verts[2];

        edges[I4(6,4,2,3, 2,1,1,0)] = verts[0];
        edges[I4(6,4,2,3, 2,1,1,1)] = verts[1];
        edges[I4(6,4,2,3, 2,1,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 2,2,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 2,2,0,1)] = verts[1];
        edges[I4(6,4,2,3, 2,2,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 2,2,1,0)] = verts[0];
        edges[I4(6,4,2,3, 2,2,1,1)] = verts[1];
        edges[I4(6,4,2,3, 2,2,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 2,3,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 2,3,0,1)] = verts[1];
        edges[I4(6,4,2,3, 2,3,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 2,3,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 2,3,1,1)] = verts[1];
        edges[I4(6,4,2,3, 2,3,1,2)] = verts[2];


        /* easT */
        edges[I4(6,4,2,3, 3,0,0,0)] = verts[0];
        edges[I4(6,4,2,3, 3,0,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,0,0,2)] = verts[2];

        edges[I4(6,4,2,3, 3,0,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 3,0,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,0,1,2)] = verts[2];


        edges[I4(6,4,2,3, 3,1,0,0)] = verts[0];
        edges[I4(6,4,2,3, 3,1,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,1,0,2)] = verts[2];

        edges[I4(6,4,2,3, 3,1,1,0)] = verts[0];
        edges[I4(6,4,2,3, 3,1,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,1,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 3,2,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 3,2,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,2,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 3,2,1,0)] = verts[0];
        edges[I4(6,4,2,3, 3,2,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,2,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 3,3,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 3,3,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,3,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 3,3,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 3,3,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 3,3,1,2)] = verts[2];


        /* froNt */
        edges[I4(6,4,2,3, 4,0,0,0)] = verts[0];
        edges[I4(6,4,2,3, 4,0,0,1)] = verts[1];
        edges[I4(6,4,2,3, 4,0,0,2)] = verts[2];

        edges[I4(6,4,2,3, 4,0,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 4,0,1,1)] = verts[1];
        edges[I4(6,4,2,3, 4,0,1,2)] = verts[2];


        edges[I4(6,4,2,3, 4,1,0,0)] = verts[0];
        edges[I4(6,4,2,3, 4,1,0,1)] = verts[1];
        edges[I4(6,4,2,3, 4,1,0,2)] = verts[2];

        edges[I4(6,4,2,3, 4,1,1,0)] = verts[0];
        edges[I4(6,4,2,3, 4,1,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 4,1,1,2)] = verts[2];


        edges[I4(6,4,2,3, 4,2,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 4,2,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 4,2,0,2)] = verts[2];

        edges[I4(6,4,2,3, 4,2,1,0)] = verts[0];
        edges[I4(6,4,2,3, 4,2,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 4,2,1,2)] = verts[2];


        edges[I4(6,4,2,3, 4,3,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 4,3,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 4,3,0,2)] = verts[2];

        edges[I4(6,4,2,3, 4,3,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 4,3,1,1)] = verts[1];
        edges[I4(6,4,2,3, 4,3,1,2)] = verts[2];


        /* bacK */
        edges[I4(6,4,2,3, 5,0,0,0)] = verts[0];
        edges[I4(6,4,2,3, 5,0,0,1)] = verts[1];
        edges[I4(6,4,2,3, 5,0,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 5,0,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 5,0,1,1)] = verts[1];
        edges[I4(6,4,2,3, 5,0,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 5,1,0,0)] = verts[0];
        edges[I4(6,4,2,3, 5,1,0,1)] = verts[1];
        edges[I4(6,4,2,3, 5,1,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 5,1,1,0)] = verts[0];
        edges[I4(6,4,2,3, 5,1,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 5,1,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 5,2,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 5,2,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 5,2,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 5,2,1,0)] = verts[0];
        edges[I4(6,4,2,3, 5,2,1,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 5,2,1,2)] = verts[2] + length;


        edges[I4(6,4,2,3, 5,3,0,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 5,3,0,1)] = verts[1] + length;
        edges[I4(6,4,2,3, 5,3,0,2)] = verts[2] + length;

        edges[I4(6,4,2,3, 5,3,1,0)] = verts[0] + length;
        edges[I4(6,4,2,3, 5,3,1,1)] = verts[1];
        edges[I4(6,4,2,3, 5,3,1,2)] = verts[2] + length;

        /* loop over faces/edges */
        for (int f = 0; f < 6; f++) {

            //printf("%ld %d\n", icell, f);

            for (int e = 0; e < 4; e++) {
                double x[3],y[3];

                if (intersection(&edges[I4(6,4,2,3,f,e,0,0)],&edges[I4(6,4,2,3,f,e,1,0)],p,n,x)) {

                    // printf("x: %ld %f %f %f\n", icell, x[0],x[1],x[2]);

                    for (int g = e; g < 4; g++) {

                        if (intersection(&edges[I4(6,4,2,3,f,g,0,0)],&edges[I4(6,4,2,3,f,g,1,0)],p,n,y)) {

                            //printf("y: %ld %f %f %f\n", icell, y[0],y[1],y[2]);

                            if (dist3(x,y) < 1e-9)
                                continue;

                            if (doedges) {
                                lines[I3(nlines[0],2,3,ecount,0,0)] = x[0];
                                lines[I3(nlines[0],2,3,ecount,0,1)] = x[1];
                                lines[I3(nlines[0],2,3,ecount,0,2)] = x[2];

                                lines[I3(nlines[0],2,3,ecount,1,0)] = y[0];
                                lines[I3(nlines[0],2,3,ecount,1,1)] = y[1];
                                lines[I3(nlines[0],2,3,ecount,1,2)] = y[2];

                                //printf("%ld %f %f %f\n", icell, x[0],x[1],x[2]);
                                //printf("%ld %f %f %f\n", icell, y[0],y[1],y[2]);
                                //printf("\n");
                            }

                            ecount++;
                        }

                    } // g
                    // printf("\n");
                }
            } // e
        } // f
    } // icell

    return ecount;
}

# if defined(P4EST)

void
morton_to_coords(
    const int dims_levels[1], const int8_t *levels,
    const int dims_morton[2], const int32_t *morton,
    const int dims_coords[2], double *coords
) {
    p4est_connectivity_t *unitcube = p4est_connectivity_new_unitsquare();

    const int nc = dims_coords[0]; // # of cells

    for (size_t icell = 0; icell < nc; icell++) {
        const double length = 1. / pow(2,levels[icell]);

        double verts[3];
        p4est_qcoord_to_vertex(unitcube, 0, 
            morton[I2(nc,2,icell,0)], morton[I2(nc,2,icell,1)], verts);

        coords[I2(nc,3,icell,0)] = verts[0];
        coords[I2(nc,3,icell,1)] = verts[1];
        coords[I2(nc,3,icell,2)] = length;
    }

    p4est_connectivity_destroy(unitcube);
}

void
cells_to_image(
    const int dims_levels[1], const int8_t *levels,
    const int dims_morton[2], const int32_t *morton,
    const int dims_cells[3], const double *cells,
    const int dims_image[2], double *const image,
    const int method, const int gridlines
) {
    p4est_connectivity_t *unitcube = p4est_connectivity_new_unitsquare();

    const int nc = dims_cells[0];
    const int nx = dims_cells[1];
    const int ny = dims_cells[2];

    const int idx = dims_image[0];
    const int idy = dims_image[1];

    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);

    for (size_t icell = 0; icell < nc; icell++) {
        const double length = 1. / pow(2,levels[icell]);

        double verts[3];
        p4est_qcoord_to_vertex(unitcube, 0, 
            morton[I2(nc,2,icell,0)], morton[I2(nc,2,icell,1)], verts);

        for (int i = 0; i < nx; i++)
            xnodes[i] = verts[0] + (i+0.5)/nx * length;

        for (int i = 0; i < ny; i++)
            ynodes[i] = verts[1] + (i+0.5)/ny * length;

        const int imgx = idx *  verts[0]           + 0.5*length/nx     + gridlines;
        const int Imgx = idx * (verts[0] + length) - 0.5*length/nx + 1 - gridlines;

        const int imgy = idy *  verts[1]           + 0.5*length/ny     + gridlines;
        const int Imgy = idy * (verts[1] + length) - 0.5*length/ny + 1 - gridlines;

        for (int i = imgx; i < Imgx; i++)
        for (int j = imgy; j < Imgy; j++)
        {
            const double x = (i+0.5)/idx;
            const double y = (j+0.5)/idy;
            
            switch(method) {
                case NEAREST:
                    image[I2(idx,idy,i,j)] = nearest2D(nx, xnodes, ny,ynodes, &cells[I3(nc,nx,ny,icell,0,0)],x,y);
                    break;

                case LINEAR:
                    image[I2(idx,idy,i,j)] = bilinear(nx, xnodes, ny,ynodes, &cells[I3(nc,nx,ny,icell,0,0)],x,y);
                    break;

                case COSINE:
                    image[I2(idx,idy,i,j)] = bicosine(nx, xnodes, ny,ynodes, &cells[I3(nc,nx,ny,icell,0,0)],x,y);
                    break;
            }
        }
    }

    free(xnodes);
    free(ynodes);
    p4est_connectivity_destroy(unitcube);
}

void
cells_to_image_3d(
    const int dims_levels[1], const int8_t *levels,
    const int dims_morton[2], const int32_t *morton,
    const int dims_cells[4], const double *cells,
    const int dims_image[3], double *const image
) {
    p8est_connectivity_t *unitcube = p8est_connectivity_new_unitcube();

    const int nc = dims_cells[0];
    const int nx = dims_cells[1];
    const int ny = dims_cells[2];
    const int nz = dims_cells[3];

    const int idx = dims_image[0];
    const int idy = dims_image[1];
    const int idz = dims_image[2];

    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);
    double *const znodes = malloc(sizeof(double) * nz);

    for (size_t icell = 0; icell < nc; icell++) {
        const double length = 1. / pow(2,levels[icell]);

        double verts[3];
        p8est_qcoord_to_vertex(unitcube, 0, 
            morton[I2(nc,3,icell,0)], morton[I2(nc,3,icell,1)], morton[I2(nc,3,icell,2)], verts);

        for (int i = 0; i < nx; i++)
            xnodes[i] = verts[0] + (i+0.5)/nx * length;

        for (int i = 0; i < ny; i++)
            ynodes[i] = verts[1] + (i+0.5)/ny * length;

        for (int i = 0; i < nz; i++)
            znodes[i] = verts[2] + (i+0.5)/nz * length;

        const int imgx = idx *  verts[0]           + 0.5*length/nx     + GRIDLINES;
        const int Imgx = idx * (verts[0] + length) - 0.5*length/nx + 1 - GRIDLINES;

        const int imgy = idy *  verts[1]           + 0.5*length/ny     + GRIDLINES;
        const int Imgy = idy * (verts[1] + length) - 0.5*length/ny + 1 - GRIDLINES;

        const int imgz = idz *  verts[2]           + 0.5*length/nz     + GRIDLINES;
        const int Imgz = idz * (verts[2] + length) - 0.5*length/nz + 1 - GRIDLINES;

        for (int i = imgx; i < Imgx; i++)
        for (int j = imgy; j < Imgy; j++)
        for (int k = imgz; k < Imgz; k++)
        {
            const double x = (i+0.5)/idx;
            const double y = (j+0.5)/idy;
            const double z = (k+0.5)/idz;

            image[I3(idx,idy,idz,i,j,k)] = nearest3D(
                nx,xnodes,ny,ynodes,nz,znodes,
                &cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);
        }
    }

    free(xnodes);
    free(ynodes);
    free(znodes);
    p8est_connectivity_destroy(unitcube);
}

# endif

void
cells_to_image_flash_ug_2d(
    const int dims_coords[2], const double *coords,
    const int dims_bsizes[2], const double *bsizes,
    const int dims_blocks[3], const double *blocks,
    const int dims_image[2], double *const image,
    const int method
){
    const int nb = dims_blocks[0]; // # of blocks
    const int nx = dims_blocks[1];
    const int ny = dims_blocks[2];

    const int ix = dims_image[0];
    const int iy = dims_image[1];

    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);

    for (size_t iblock = 0; iblock < nb; iblock++) {
        const double xlen = bsizes[I2(nb,2,iblock,0)];
        const double ylen = bsizes[I2(nb,2,iblock,1)];

        const double xmid = coords[I2(nb,2,iblock,0)];
        const double ymid = coords[I2(nb,2,iblock,1)];

        for (int i = 0; i < nx; i++)
            xnodes[i] = xmid - 0.5*xlen + (i+0.5)*xlen/nx;

        for (int i = 0; i < ny; i++)
            ynodes[i] = ymid - 0.5*ylen + (i+0.5)*ylen/ny;

        // map block vertices onto image space
        const int imgx = ix * (xmid - 0.5*xlen) + 0.5*xlen/nx;
        const int Imgx = ix * (xmid + 0.5*xlen) + 0.5*xlen/nx;

        const int imgy = iy * (ymid - 0.5*ylen) + 0.5*ylen/ny;
        const int Imgy = iy * (ymid + 0.5*ylen) + 0.5*ylen/ny;

        for (int i = imgx; i < Imgx; i++)
        for (int j = imgy; j < Imgy; j++)
        {
            const double x = (i+0.5)/ix;
            const double y = (j+0.5)/iy;
            
            switch(method) {
                case NEAREST:
                    image[I2(ix,iy,i,j)] = nearest2D(nx, xnodes, ny,ynodes, &blocks[I3(nb,nx,ny,iblock,0,0)],x,y);
                    break;

                case LINEAR:
                    image[I2(ix,iy,i,j)] = bilinear(nx, xnodes, ny,ynodes, &blocks[I3(nb,nx,ny,iblock,0,0)],x,y);
                    break;

                case COSINE:
                    image[I2(ix,iy,i,j)] = bicosine(nx, xnodes, ny,ynodes, &blocks[I3(nb,nx,ny,iblock,0,0)],x,y);
                    break;
            }
        }
    }

    free(xnodes);
    free(ynodes);
}

void
cells_to_image_titanic_patch_2d(
    const int pshape[4], const double *patch,
    const int ishape[2], double *const image,
    const int method
){
    const double DX = 1.0;
    const double DY = 1.0;

    const int Nx = pshape[0];
    const int Ny = pshape[1];

    const int nx = pshape[2];
    const int ny = pshape[3];

    const int Ix = ishape[0];
    const int Iy = ishape[1];

    //printf("%d %d | %d %d | %d %d\n",Nx,Ny,nx,ny,Ix,Iy);

    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);

    const double dx = DX/((double) Nx);
    const double dy = DY/((double) Ny);

    const double Dx = DX/((double) Ix);
    const double Dy = DY/((double) Iy);

    for (int px = 0; px < Nx; px++) {
        const double x0 = dx*px;

        for (int i = 0; i < nx; i++)
            xnodes[i] = x0 + (i+0.5)*dx/nx;

        // map patch space into image space
        const int I0 = Ix *  x0;
        const int I1 = Ix * (x0 + dx);

        for (int py = 0; py < Ny; py++) {
            const double y0 = dy*py;

            for (int i = 0; i < ny; i++)
                ynodes[i] = y0 + (i+0.5)*dy/ny;

            const int J0 = Iy *  y0;
            const int J1 = Iy * (y0 + dy);

            //printf("%d %d | %f %f | %d %d | %d %d\n",px,py,x0,y0,I0,I1,J0,J1);
            //printf("%d %d ",px,py);
            //printf("| %f %f | %f %f ",xnodes[0],xnodes[nx-1],Dx*(I0+0.5),Dx*(I1+0.5));
            //printf("| %f %f | %f %f\n",ynodes[0],ynodes[ny-1],Dy*(J0+0.5),Dy*(J1+0.5));

            //printf("%f %f | %d %d | %d %d | %d %d\n",x0,y0,px,py,I0,I1,J0,J1);

            for (int i = I0; i < I1; i++)
            for (int j = J0; j < J1; j++)
            {
                const double x = Dx*(i+0.5);
                const double y = Dy*(j+0.5);
                
                //printf("%d %d | %f %f | %f %f\n",i,j,x,y,xnodes[0],xnodes[nx-1]);

                switch(method) {
                    case NEAREST:
                        image[I2(Ix,Iy,i,j)] = nearest2D(nx, xnodes, ny,ynodes, &patch[I4(Nx,Ny,nx,ny,px,py,0,0)],x,y);
                        break;

                    case LINEAR:
                        image[I2(Ix,Iy,i,j)] = bilinear(nx, xnodes, ny,ynodes, &patch[I4(Nx,Ny,nx,ny,px,py,0,0)],x,y);
                        break;
                }
            }
        }
    }

    free(xnodes);
    free(ynodes);
}

/* ========================================================================= */
# if defined(P4EST)

void
cells_to_plane_3d(
    const int dlevels[1], const int8_t *levels,
    const int dmorton[2], const int32_t *morton,
    const int dcells[4], const double *cells,
    const int dimage[2], double *const image,
    const double p[3], const double u[3], const double v[3],
    const int method
) {
    p8est_connectivity_t *unitcube = p8est_connectivity_new_unitcube();

    const int nc = dcells[0];
    const int nx = dcells[1];
    const int ny = dcells[2];
    const int nz = dcells[3];

    // physical locations of nodes within cube
    double *const xnodes = malloc(sizeof(double) * nx);
    double *const ynodes = malloc(sizeof(double) * ny);
    double *const znodes = malloc(sizeof(double) * nz);

    const double uu = dot3(u,u);
    const double vv = dot3(v,v);

    const double su = sqrt(uu);
    const double sv = sqrt(vv);

    double normal[3];
    cross3(u,v,normal);

    for (size_t icell = 0; icell < nc; icell++) {

        double verts[3];
        double center[3];
        double pclose[3];
        double vshift[3];

        const double length = 1. / pow(2,levels[icell]);
        const double half = 0.5*length;

        p8est_qcoord_to_vertex(unitcube, 0, 
            morton[I2(nc,3,icell,0)], morton[I2(nc,3,icell,1)], morton[I2(nc,3,icell,2)], verts);

        /* ----------------------------------------------------------------- */

        center[0] = verts[0] + half;
        center[1] = verts[1] + half;
        center[2] = verts[2] + half;

        pclose3(p,normal,center,pclose);

        if (dist3(center,pclose) > length)
            continue;
            
        /* ----------------------------------------------------------------- */

        for (int i = 0; i < nx; i++)
            xnodes[i] = verts[0] + (i+0.5)/nx * length;

        for (int i = 0; i < ny; i++)
            ynodes[i] = verts[1] + (i+0.5)/ny * length;

        for (int i = 0; i < nz; i++)
            znodes[i] = verts[2] + (i+0.5)/nz * length;

        /* ----------------------------------------------------------------- */

# if 0
        // vshift[0] = verts[0]-p[0];
        // vshift[1] = verts[1]-p[1];
        // vshift[2] = verts[2]-p[2];

        vshift[0] = pclose[0] - 2*length*(u[0]/su + v[0]/sv);
        vshift[1] = pclose[1] - 2*length*(u[1]/su + v[1]/sv);
        vshift[2] = pclose[2] - 2*length*(u[2]/su + v[2]/sv);

        const int ix = imax(0,dimage[0] * dot3(vshift,u)/uu - 0.5);
        const int iy = imax(0,dimage[1] * dot3(vshift,v)/vv - 0.5);

        /* ----------------------------------------------------------------- */

        // vshift[0] = verts[0]-p[0] + length;
        // vshift[1] = verts[1]-p[1] + length;
        // vshift[2] = verts[2]-p[2] + length;

        vshift[0] = pclose[0] + 2*length*(u[0]/su + v[0]/sv);
        vshift[1] = pclose[1] + 2*length*(u[1]/su + v[1]/sv);
        vshift[2] = pclose[2] + 2*length*(u[2]/su + v[2]/sv);

        const int Ix = imin(dimage[0]-1,dimage[0] * dot3(vshift,u)/uu - 0.5);
        const int Iy = imin(dimage[0]-1,dimage[1] * dot3(vshift,v)/vv - 0.5);

        //printf("%d %d %d %d\n", ix, Ix, iy, Iy);
# else
        const int ix = 0;
        const int iy = 0;

        const int Ix = dimage[0]-1;
        const int Iy = dimage[1]-1;
# endif

        /* ----------------------------------------------------------------- */

        for (int i = ix; i <= Ix; i++)
        for (int j = iy; j <= Iy; j++)
        {
            const double x = p[0] + (i+0.5)/dimage[0]*u[0] + (j+0.5)/dimage[1]*v[0];
            const double y = p[1] + (i+0.5)/dimage[0]*u[1] + (j+0.5)/dimage[1]*v[1];
            const double z = p[2] + (i+0.5)/dimage[0]*u[2] + (j+0.5)/dimage[1]*v[2];

            if (!(fabs(center[0]-x) <= half
               && fabs(center[1]-y) <= half
               && fabs(center[2]-z) <= half)) {
                continue;
            }

            switch(method) {
                case NEAREST:
                    image[I2(dimage[0],dimage[1],i,j)]
                        = nearest3D(nx,xnodes,ny,ynodes,nz,znodes,&cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);
                    break;

                case LINEAR:
                    image[I2(dimage[0],dimage[1],i,j)]
                        = trilinear(nx,xnodes,ny,ynodes,nz,znodes,&cells[I4(nc,nx,ny,nz,icell,0,0,0)],x,y,z);
                    break;
            }
        }
    }

    free(xnodes);
    free(ynodes);
    free(znodes);
    p8est_connectivity_destroy(unitcube);
}

# endif

# if defined(P4EST)

int
plane_morton_to_coords(
    const int dlevels[1], const int8_t *levels,
    const int dmorton[2], const int32_t *morton,
    const double p[3], const double u[3], const double v[3],
    const int dedges[3], double *edges, const int doedges
) {
    p8est_connectivity_t *unitcube = p8est_connectivity_new_unitcube();

    const int nc = dlevels[0]; // # of cells

    const double uu = dot3(u,u);
    const double vv = dot3(v,v);

    const double lu = sqrt(uu);
    const double lv = sqrt(vv);

    double n[3];
    cross3(u,v,n);

    int ecount = 0;

    for (size_t icell = 0; icell < nc; icell++) {
        double verts[3];
        double center[3];
        double pclose[3];
        double vshift[3];
        double temp;

        const double length = 1. / pow(2,levels[icell]);
        const double half = 0.5*length;

        double edges[6*4*2*3];

        p8est_qcoord_to_vertex(unitcube, 0, 
            morton[i2(nc,3,icell,0)], morton[i2(nc,3,icell,1)], morton[i2(nc,3,icell,2)], verts);

        /* ----------------------------------------------------------------- */

        center[0] = verts[0] + half;
        center[1] = verts[1] + half;
        center[2] = verts[2] + half;

        pclose3(p,n,center,pclose);

        if (!(fabs(center[0]-pclose[0]) <= half
           && fabs(center[1]-pclose[1]) <= half
           && fabs(center[2]-pclose[2]) <= half)) {
            continue;
        }

        vshift[0] = pclose[0] - p[0];
        vshift[1] = pclose[1] - p[1];
        vshift[2] = pclose[2] - p[2];

        //temp = dot3(vshift,u);
        //if (!(0 <= temp && temp <= uu))
        //    continue;
 
        //temp = dot3(vshift,v);
        //if (!(0 <= temp && temp <= vv))
        //    continue;
 
        /* north */
        edges[i4(6,4,2,3, 0,0,0,0)] = verts[0];
        edges[i4(6,4,2,3, 0,0,0,1)] = verts[1];
        edges[i4(6,4,2,3, 0,0,0,2)] = verts[2];

        edges[i4(6,4,2,3, 0,0,1,0)] = verts[0];
        edges[i4(6,4,2,3, 0,0,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 0,0,1,2)] = verts[2];


        edges[i4(6,4,2,3, 0,1,0,0)] = verts[0];
        edges[i4(6,4,2,3, 0,1,0,1)] = verts[1];
        edges[i4(6,4,2,3, 0,1,0,2)] = verts[2];

        edges[i4(6,4,2,3, 0,1,1,0)] = verts[0];
        edges[i4(6,4,2,3, 0,1,1,1)] = verts[1];
        edges[i4(6,4,2,3, 0,1,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 0,2,0,0)] = verts[0];
        edges[i4(6,4,2,3, 0,2,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 0,2,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 0,2,1,0)] = verts[0];
        edges[i4(6,4,2,3, 0,2,1,1)] = verts[1];
        edges[i4(6,4,2,3, 0,2,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 0,3,0,0)] = verts[0];
        edges[i4(6,4,2,3, 0,3,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 0,3,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 0,3,1,0)] = verts[0];
        edges[i4(6,4,2,3, 0,3,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 0,3,1,2)] = verts[2];


        /* south */
        edges[i4(6,4,2,3, 1,0,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,0,0,1)] = verts[1];
        edges[i4(6,4,2,3, 1,0,0,2)] = verts[2];

        edges[i4(6,4,2,3, 1,0,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,0,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 1,0,1,2)] = verts[2];


        edges[i4(6,4,2,3, 1,1,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,1,0,1)] = verts[1];
        edges[i4(6,4,2,3, 1,1,0,2)] = verts[2];

        edges[i4(6,4,2,3, 1,1,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,1,1,1)] = verts[1];
        edges[i4(6,4,2,3, 1,1,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 1,2,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,2,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 1,2,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 1,2,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,2,1,1)] = verts[1];
        edges[i4(6,4,2,3, 1,2,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 1,3,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,3,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 1,3,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 1,3,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 1,3,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 1,3,1,2)] = verts[2];



        /* west */
        edges[i4(6,4,2,3, 2,0,0,0)] = verts[0];
        edges[i4(6,4,2,3, 2,0,0,1)] = verts[1];
        edges[i4(6,4,2,3, 2,0,0,2)] = verts[2];

        edges[i4(6,4,2,3, 2,0,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 2,0,1,1)] = verts[1];
        edges[i4(6,4,2,3, 2,0,1,2)] = verts[2];


        edges[i4(6,4,2,3, 2,1,0,0)] = verts[0];
        edges[i4(6,4,2,3, 2,1,0,1)] = verts[1];
        edges[i4(6,4,2,3, 2,1,0,2)] = verts[2];

        edges[i4(6,4,2,3, 2,1,1,0)] = verts[0];
        edges[i4(6,4,2,3, 2,1,1,1)] = verts[1];
        edges[i4(6,4,2,3, 2,1,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 2,2,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 2,2,0,1)] = verts[1];
        edges[i4(6,4,2,3, 2,2,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 2,2,1,0)] = verts[0];
        edges[i4(6,4,2,3, 2,2,1,1)] = verts[1];
        edges[i4(6,4,2,3, 2,2,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 2,3,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 2,3,0,1)] = verts[1];
        edges[i4(6,4,2,3, 2,3,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 2,3,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 2,3,1,1)] = verts[1];
        edges[i4(6,4,2,3, 2,3,1,2)] = verts[2];


        /* east */
        edges[i4(6,4,2,3, 3,0,0,0)] = verts[0];
        edges[i4(6,4,2,3, 3,0,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,0,0,2)] = verts[2];

        edges[i4(6,4,2,3, 3,0,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 3,0,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,0,1,2)] = verts[2];


        edges[i4(6,4,2,3, 3,1,0,0)] = verts[0];
        edges[i4(6,4,2,3, 3,1,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,1,0,2)] = verts[2];

        edges[i4(6,4,2,3, 3,1,1,0)] = verts[0];
        edges[i4(6,4,2,3, 3,1,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,1,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 3,2,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 3,2,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,2,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 3,2,1,0)] = verts[0];
        edges[i4(6,4,2,3, 3,2,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,2,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 3,3,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 3,3,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,3,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 3,3,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 3,3,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 3,3,1,2)] = verts[2];


        /* front */
        edges[i4(6,4,2,3, 4,0,0,0)] = verts[0];
        edges[i4(6,4,2,3, 4,0,0,1)] = verts[1];
        edges[i4(6,4,2,3, 4,0,0,2)] = verts[2];

        edges[i4(6,4,2,3, 4,0,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 4,0,1,1)] = verts[1];
        edges[i4(6,4,2,3, 4,0,1,2)] = verts[2];


        edges[i4(6,4,2,3, 4,1,0,0)] = verts[0];
        edges[i4(6,4,2,3, 4,1,0,1)] = verts[1];
        edges[i4(6,4,2,3, 4,1,0,2)] = verts[2];

        edges[i4(6,4,2,3, 4,1,1,0)] = verts[0];
        edges[i4(6,4,2,3, 4,1,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 4,1,1,2)] = verts[2];


        edges[i4(6,4,2,3, 4,2,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 4,2,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 4,2,0,2)] = verts[2];

        edges[i4(6,4,2,3, 4,2,1,0)] = verts[0];
        edges[i4(6,4,2,3, 4,2,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 4,2,1,2)] = verts[2];


        edges[i4(6,4,2,3, 4,3,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 4,3,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 4,3,0,2)] = verts[2];

        edges[i4(6,4,2,3, 4,3,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 4,3,1,1)] = verts[1];
        edges[i4(6,4,2,3, 4,3,1,2)] = verts[2];


        /* back */
        edges[i4(6,4,2,3, 5,0,0,0)] = verts[0];
        edges[i4(6,4,2,3, 5,0,0,1)] = verts[1];
        edges[i4(6,4,2,3, 5,0,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 5,0,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 5,0,1,1)] = verts[1];
        edges[i4(6,4,2,3, 5,0,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 5,1,0,0)] = verts[0];
        edges[i4(6,4,2,3, 5,1,0,1)] = verts[1];
        edges[i4(6,4,2,3, 5,1,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 5,1,1,0)] = verts[0];
        edges[i4(6,4,2,3, 5,1,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 5,1,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 5,2,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 5,2,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 5,2,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 5,2,1,0)] = verts[0];
        edges[i4(6,4,2,3, 5,2,1,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 5,2,1,2)] = verts[2] + length;


        edges[i4(6,4,2,3, 5,3,0,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 5,3,0,1)] = verts[1] + length;
        edges[i4(6,4,2,3, 5,3,0,2)] = verts[2] + length;

        edges[i4(6,4,2,3, 5,3,1,0)] = verts[0] + length;
        edges[i4(6,4,2,3, 5,3,1,1)] = verts[1];
        edges[i4(6,4,2,3, 5,3,1,2)] = verts[2] + length;

        /* loop over faces/edges */
        for (int f = 0; f < 6; f++) {
            for (int e = 0; e < 4; e++) {
                for (int e = e; e < 4; e++) {
                    double x[3],y[3];

                    if(intersection(&edges[i4(6,4,2,3,f,e,0,0)],&edges[i4(6,4,2,3,f,e,1,0)],p,n,x) &&
                       intersection(&edges[i4(6,4,2,3,f,e,0,0)],&edges[i4(6,4,2,3,f,e,1,0)],p,n,y)) {

                        if (dist3(x,y) < 1e-9)
                            continue;

                        if (doedges) {
                            edges[i3(dedges[0],2,3,ecount,0,0)] = x[0];
                            edges[i3(dedges[0],2,3,ecount,0,1)] = x[1];
                            edges[i3(dedges[0],2,3,ecount,0,2)] = x[2];

                            edges[i3(dedges[0],2,3,ecount,1,0)] = y[0];
                            edges[i3(dedges[0],2,3,ecount,1,1)] = y[1];
                            edges[i3(dedges[0],2,3,ecount,1,2)] = y[2];

                            //printf("%ld %f %f %f\n", icell, x[0],x[1],x[2]);
                            //printf("%ld %f %f %f\n", icell, y[0],y[1],y[2]);
                            //printf("\n");
                        }

                        ecount++;
                    }
                }
            }
        }
    }

    p8est_connectivity_destroy(unitcube);

    return ecount;
}

# endif
