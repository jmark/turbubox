# include <stdio.h>
# include <math.h>

//  X[Nx*Ny*Nz] -> input array
// cs[nsamples] -> counts
// rs[nsamples] -> radii
// ts[nsamples] -> totals

void shell_avg_2d(
    const double *X, const int Nx, const int Ny,
    double *cs, double *rs, double *ts, const int nsamples, const int mult_with_r
)
{
    // cell half-width offsets
    const double cxofs = .5/Nx;
    const double cyofs = .5/Ny;

    const double I_max = Nx/2. + cxofs;
    const double J_max = Ny/2. + cyofs;
    const double r_max = sqrt(I_max*I_max + J_max*J_max);

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    {
        const double I = i - Nx/2. + cxofs;
        const double J = j - Ny/2. + cyofs;
        const double r = sqrt(I*I + J*J);
        const    int R = r/r_max * nsamples;

        cs[R] += 1;
        rs[R] += r;

        if (mult_with_r) {
            ts[R] += r*X[(i * Ny) + j];
        } else {
            ts[R] +=   X[(i * Ny) + j];
        }
    }
}

void shell_avg_2d_corner(
    const double *X, const int Nx, const int Ny,
    double *cs, double *rs, double *ts, const int nsamples, const int mult_with_r
)
{
    // cell half-width offsets
    const double cxofs = .5/Nx;
    const double cyofs = .5/Ny;

    const double I_max = Nx - cxofs;
    const double J_max = Ny - cyofs;
    const double r_max = sqrt(I_max*I_max + J_max*J_max);

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    {
        const double I = i + cxofs;
        const double J = j + cyofs;
        const double r = sqrt(I*I + J*J);
        const    int R = r/r_max * nsamples;

        cs[R] += 1;
        rs[R] += r;

        if (mult_with_r) {
            ts[R] += r*X[(i * Ny) + j];
        } else {
            ts[R] +=   X[(i * Ny) + j];
        }
    }
}

void shell_avg_3d(
    const double *X, const int Nx, const int Ny, const int Nz,
    double *cs, double *rs, double *ts, const int nsamples, const int mult_with_rsquare
)
{
    // cell half-width offsets
    const double cxofs = .5/Nx;
    const double cyofs = .5/Ny;
    const double czofs = .5/Nz;

    const double I_max = Nx/2. + cxofs;
    const double J_max = Ny/2. + cyofs;
    const double K_max = Nz/2. + czofs;
    const double r_max = sqrt(I_max*I_max + J_max*J_max + K_max*K_max);

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        const double I = i - Nx/2. + cxofs;
        const double J = j - Ny/2. + cyofs;
        const double K = k - Nz/2. + czofs;
        const double r = sqrt(I*I + J*J + K*K);
        const    int R = r/r_max * nsamples;

        cs[R] += 1;
        rs[R] += r;

        if (mult_with_rsquare) {
            ts[R] += r*r*X[((i * Ny) + j) * Nz + k];
        } else {
            ts[R] += X[((i * Ny) + j) * Nz + k];
        }
    }
}

void shell_max_3d(
    const double *X, const int Nx, const int Ny, const int Nz,
    double *cs, double *rs, double *ts, const int nsamples
)
{
    // cell half-width offsets
    const double cxofs = .5/Nx;
    const double cyofs = .5/Ny;
    const double czofs = .5/Nz;

    const double I_max = Nx/2. + cxofs;
    const double J_max = Ny/2. + cyofs;
    const double K_max = Nz/2. + czofs;
    const double r_max = sqrt(I_max*I_max + J_max*J_max + K_max*K_max);

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        const double I = i - Nx/2. + cxofs;
        const double J = j - Ny/2. + cyofs;
        const double K = k - Nz/2. + czofs;
        const double r = sqrt(I*I + J*J + K*K);
        const    int R = r/r_max * nsamples;

        cs[R] += 1;
        rs[R] += r;

        const double old = ts[R];
        const double new = X[((i * Ny) + j) * Nz + k];

        ts[R] = new > old ? new : old;
    }
}
