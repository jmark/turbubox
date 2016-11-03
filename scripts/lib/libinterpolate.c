# include <stdio.h>
# include <math.h>

double LagrangePolynomial(const double *xs, const int xslen, const int j, const double x) {
    double acc = 1;
    
    for (int i = 0; i < xslen; ++i) {
        if (i==j) continue;
        acc *= (x-xs[i])/(xs[j]-xs[i]);
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
    const int *Ks, 

    const int *IOs, 
    const int *JOs, 
    const int *KOs, const int IJKslen,

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
