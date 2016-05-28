#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <cmath>

# include "H5File.hpp"
# include "H5DataSet.hpp"
# include "H5DataSpace.hpp"
# include "H5DataType.hpp"

# include "Array.h"
# include "fftw++.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace std;
using namespace utils;
using namespace Array;
using namespace fftwpp;

typedef unsigned int uint;
typedef std::vector<double> dvec;

double norm3d (const double x, const double y, const double z)
{
    return std::sqrt( x*x + y*y + z*z );
}

double norm3d (const dvec & v)
{
    return std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

Array::array3<double> &read_ds
(
    const HighFive::File &file,
    const std::string &dname,
    const std::vector<uint> &dims
)
{
    auto dataset = file.getDataSet(dname);

    std::vector<double> tmp;
    dataset.read(tmp);

    auto *arr = new Array::array3<double>(dims[0],dims[1],dims[2]);
    *arr = tmp.data();

    return *arr;
}

void read_meta
(
    const HighFive::File &file,
    const std::string &dname,
    std::vector<uint> &vec
)
{
    auto dataset = file.getDataSet(dname);
    dataset.read(vec);
}

int main(int argc, char * argv[])
{
    // ---------------------------------------------------------------------- //
    // Read in Commandline args

    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " filename" << endl;
        return 1;
    }

    int argPtr = 1;
    const std::string file_name = argv[argPtr++];

    // ---------------------------------------------------------------------- //
    // Initialize arrays for FFTW++

    // Open input file
    HighFive::File h5file(file_name, HighFive::File::ReadOnly);

    // Get Meta Data
    std::vector<uint> dims(3);
    read_meta(h5file,"DIMS",dims);

    cerr 
        << dims[0] << "\t"
        << dims[1] << "\t"
        << dims[2] << "\t"
        << endl;

    auto &dens = read_ds(h5file,"dens",dims);

    cerr 
        << dens(0,0,0) << "\t"
        << dens(0,30,1) << "\t"
        << dens(63,63,63) << "\t"
        << endl;

    // fftw::maxthreads = get_max_threads();
    fftw::maxthreads = 1;

    unsigned int Nzp = dims[0]/2+1;
    size_t align = sizeof(Complex);

    // in
    //array3<double> dens_r(Nx,Ny,Nz,align);
    //array3<double> ergy_r(Nx,Ny,Nz,align);

    // out
    array3<Complex> dens_c(dims[0],dims[1],Nzp,align);

    // process
    rcfft3d dens_forward(dims[0],dims[1],dims[2],dens,dens_c);

    // ---------------------------------------------------------------------- //
    // Compute Fourier Spectra

    dens_forward.fft0Normalized(dens,dens_c);

    //  vector<double> dens_avg(Nx);
    //  vector<double> ergy_avg(Nx);
    //  vector<uint> avg_n(Nx);


    //  // ---------------------------------------------------------------------- //
    //  // Calculate spherical shell average

    //  for (uint i = 0 ; i < Nx ; i++)
    //  {
    //      for (uint j = 0 ; j < Ny ; j++)
    //      {
    //          for (uint k = 0 ; k < Nzp ; k++)
    //          {
    //              int p = i - Nx/2;
    //              int q = j - Ny/2;
    //              int r = k;

    //              double K = norm3d(p,q,r);
    //              dens_avg[round(K)] += K*K * abs(dens_c(i,j,k));
    //              ergy_avg[round(K)] += K*K * abs(ergy_c(i,j,k));
    //              avg_n[round(K)]++;
    //          }
    //      }
    //  }

    //  // ---------------------------------------------------------------------- //
    //  // Calculate power spectrum

    //  for (uint K = 0 ; K < Nx ; K++)
    //  {
    //      if ( avg_n[K] > 0 )
    //          cout
    //              << K << "\t"
    //              << dens_avg[K]/avg_n[K] << "\t"
    //              << ergy_avg[K]/avg_n[K] << "\t"
    //              << endl;
    //      else
    //          cout << K << 0 << 0 << endl;
    //  }

    //  return 0;

    //  // ---------------------------------------------------------------------- //
    //  // Spit out results

    //  for (uint i = 0 ; i < Nx ; i++)
    //  {
    //      double x = (i + 0.5)*dx;

    //      for (uint j = 0 ; j < Ny ; j++)
    //      {
    //          double y = (j + 0.5)*dy;

    //          for (uint k = 0 ; k < Nzp ; k++)
    //          {
    //              double z = (k + 0.5)*dz;

    //              cout
    //                  << i << "\t"
    //                  << j << "\t"
    //                  << k << "\t"

    //                  << x << "\t"
    //                  << y << "\t"
    //                  << z << "\t"

    //                  << abs( dens_c(i,j,k) ) << "\t"
    //                  << arg( dens_c(i,j,k) ) << "\t"

    //                  << abs( ergy_c(i,j,k) ) << "\t"
    //                  << arg( ergy_c(i,j,k) ) << "\t"

    //                  << endl;
    //          }
    //      }
    //      cout << endl;
    //  }

    //  return 0;
}
