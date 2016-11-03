#include <stdio.h>

int flash_to_flexi(
    const int *Is, 
    const int *Js, 
    const int *Ks, const int IJKslen,

    const double *xs, const int xslen,
    const double *ys, const int yslen,
    const double *zs, const int zslen,
    const double *fss, const int fsswidth, const int fssheight, const int fssdepth, const int fslen,

    const double *Xs, const double *Ys, const double *Zs, double *Fss, const int Fslen);

# define N_POLY     3
# define N_ITPL     5

# define xs_LEN      (N_ITPL + 1)
# define ys_LEN      (N_ITPL + 1) 
# define zs_LEN      (N_ITPL + 1) 

# define Xs_LEN      (N_POLY + 1)
# define Ys_LEN      (N_POLY + 1) 
# define Zs_LEN      (N_POLY + 1) 

# define fss_WIDTH   (128+2)
# define fss_HEIGHT  (128+2)
# define fss_DEPTH   (128+2)

# define fss_LEN     (fss_WIDTH * fss_HEIGHT * fss_DEPTH)
# define IJKs_LEN    (fss_LEN /xs_LEN/ys_LEN/zs_LEN)

int main(void)
{
    
    int Is[IJKs_LEN]; 
    int Js[IJKs_LEN]; 
    int Ks[IJKs_LEN]; 

    double xs[XSLEN];
    double ys[YSLEN];
    double zs[ZSLEN];
    
    double fss[IJKs_LEN * FSLEN];
}
