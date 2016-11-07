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
    const int Xslen, const double *Xs,       double *Fss,
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
            Fs[((I * Ny) + J) * Nz + K] = F;
            //printf("%d %d %d -> %f\n", I,J,K, F);
        }
    }

    free(Ls);
}
