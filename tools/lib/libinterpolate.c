# include <stdio.h>
# include <math.h>
# include <malloc.h>

double LagrangePolynomial(const double *xs, const int xslen, const int j, const double X) {
    double acc = 1;
    
    for (int i = 0; i < xslen; ++i) {
        if (i==j) continue;
        acc *= (X-xs[i])/(xs[j]-xs[i]);
    }

    return acc;
}

int lagrange_interpolate_3d(
    const double *xs, const int xslen,
    const double *ys, const int yslen,
    const double *zs, const int zslen,
    const double *fs, const int fslen,
    const double *Xs, const double *Ys, const double *Zs, double *Fs, const int Fslen)
{
    for (int m = 0; m < Fslen; m++)
    {
        double X = Xs[m];
        double Y = Ys[m];
        double Z = Zs[m];
        double F = 0;

        for (int i = 0; i < xslen; i++)
        {
            const double Lx = LagrangePolynomial(xs, xslen, i, X);
            for (int j = 0; j < yslen; j++)
            {
                const double Ly = LagrangePolynomial(ys, yslen, j, Y);
                for (int k = 0; k < zslen; k++)
                {
                    const double Lz = LagrangePolynomial(zs, zslen, k, Z);
                    const double f = fs[((i * yslen) + j) * zslen + k];
                    F += f*Lx*Ly*Lz;
                }
            }
        }
        Fs[m] = F;
    }
}


int flash_to_flexi(
    const int *Is, 
    const int *Js, 
    const int *Ks, const int IJKslen,

    const double *xs, const int xslen,
    const double *ys, const int yslen,
    const double *zs, const int zslen,
    const double *fss, const int fsswidth, const int fssheight, const int fssdepth, const int fslen,

    const double *Xs, const double *Ys, const double *Zs, double *Fss, const int Fslen)
{
    for (int elemid = 0; elemid < IJKslen; elemid++)
    {
        const int i = Is[elemid];
        const int j = Js[elemid];
        const int k = Ks[elemid];

        const double *fs = fss + ((i * fssheight) + j) * fssdepth + k;
              double *Fs = Fss + elemid * Fslen;

        lagrange_interpolate_3d(xs,xslen, ys,yslen, zs,zslen, fs,fslen, Xs,Ys,Zs,Fs,Fslen);
    }
}

int lagrange_interpolate_3d_RG(
    const int xslen, const int Xslen,
    const double *Ls, const double *fs, double *Fs
)
{
    for (int I = 0; I < Xslen; I++)
    for (int J = 0; J < Xslen; J++)
    for (int K = 0; K < Xslen; K++)
    {
        double F = 0;
        for (int i = 0; i < xslen; i++)
        for (int j = 0; j < xslen; j++)
        for (int k = 0; k < xslen; k++)
            F += fs[((i * xslen) + j) * xslen + k] * Ls[I*Xslen + i]*Ls[J*Xslen + j]*Ls[K*Xslen + k];

        Fs[((I * Xslen) + J) * Xslen + K] = F;
    }
}

int flash_to_flexi_RG(
    const double *xs, const int xslen, const double *fss,
    const double *Xs, const int Xslen,       double *Fss,
    const int *offsets, const int offsetslen
)
{
    double *Ls = malloc(sizeof(double) * Xslen * xslen);

    for (int I = 0; I < Xslen; I++)
    for (int i = 0; i < xslen; i++)
        Ls[I*Xslen + i] = LagrangePolynomial(xs, xslen, i, Xs[I]);

    const int Offset = Xslen * Xslen * Xslen;

    for (int elemid = 0; elemid < offsetslen; elemid++)
    {
        const double *fs = fss + offsets[elemid]; // non-consecutive
              double *Fs = Fss + Offset *elemid;  // consecutive

        lagrange_interpolate_3d_RG(xslen,Xslen, Ls,fs,Fs);
    }
}
