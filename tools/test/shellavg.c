# define XLEN 100
# define YLEN 100
# define ZLEN 100
# define SLEN 20

void shell_avg_3d(
    const double *X, const int Nx, const int Ny, const int Nz,
    double *rs, double *cs, double *ts, const int nsamples
);

int main(void)
{
    double X[XLEN * YLEN * ZLEN];
    double rs[SLEN];
    double cs[SLEN];
    double ts[SLEN];
    
    shell_avg_3d(
        X, XLEN, YLEN, ZLEN,
        rs, cs, ts, SLEN
    );
}
