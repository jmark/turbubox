#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "quickflash_file_datafile.hpp"
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_file_meshinfo.hpp"

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
    const std::string file_name     = argv[argPtr++];

    // ---------------------------------------------------------------------- //
    // Open the file
    QuickFlash::File::DataFile dfile(file_name);

    // ---------------------------------------------------------------------- //
    // Get mesh info
    const QuickFlash::File::MeshInfo meshinfo = dfile.get_mesh_info();

    // ---------------------------------------------------------------------- //
    // Get physical cell volume (uniform grid) in all directions
    const dvec min = meshinfo.get_volume_minbounds();
    const dvec max = meshinfo.get_volume_maxbounds();

    // TODO: find "smallest" cell in AMR grid
    const QuickFlash::Block::BlockInfo binfo = meshinfo.get_block_info(min);
    const dvec cvl = binfo.get_cell_width();

    const double dx = cvl[0];
    const double dy = cvl[1];
    const double dz = cvl[2];

    // ---------------------------------------------------------------------- //
    // Get number of cells (uniform grid) in all directions
    const dvec vol = meshinfo.get_volume_width();

    const uint Nx = std::round(vol[0]/cvl[0]);
    const uint Ny = std::round(vol[1]/cvl[1]);
    const uint Nz = std::round(vol[2]/cvl[2]);

    // ---------------------------------------------------------------------- //
    // DIAGNOSTICS

    // const uint nblocks = meshinfo.get_num_blocks();
    // const uint dims    = meshinfo.get_dims();

    // cerr << "Dims [ "       << dims         << " ]" << endl;
    // cerr << "Num blocks [ " << nblocks   << " ]" << endl;
    // cerr << "Num block cells [ " << meshinfo.get_num_block_cells() << " ]" << endl;

    // const std::vector<uint> block_dims        = meshinfo.get_block_dims();
    // const std::vector<uint> base_block_dims   = meshinfo.get_base_block_dims();

    // cerr << "Block dimensions [ ";
    // for (uint i = 0; i < dims; i++)
    //     cerr << block_dims[i] << " ";
    // cerr << "]" << endl;

    // cerr << "Base block dimensions [ ";
    // for (uint i = 0; i < dims; i++)
    //     cerr << base_block_dims[i] << " ";
    // cerr << "]" << endl;

    // cerr << "Min  Bounds [ " << min[0] << " " << min[1] << " " << min[2] << " ]" << endl;
    // cerr << "Max  Bounds [ " << max[0] << " " << max[1] << " " << max[2] << " ]" << endl;
    // cerr << "Vol  Widths [ " << vol[0] << " " << vol[1] << " " << vol[2] << " ]" << endl;
    // cerr << "Cell Widths [ " << cvl[0] << " " << cvl[1] << " " << cvl[2] << " ]" << endl;

    // ---------------------------------------------------------------------- //
    // Open datasets
    const QuickFlash::File::Dataset & ddens = dfile.get_dataset("dens");
    const QuickFlash::File::Dataset & dvelx = dfile.get_dataset("velx");
    const QuickFlash::File::Dataset & dvely = dfile.get_dataset("vely");
    const QuickFlash::File::Dataset & dvelz = dfile.get_dataset("velz");

    // ---------------------------------------------------------------------- //
    // Initialize arrays for FFTW++

    // fftw::maxthreads = get_max_threads();
    fftw::maxthreads = 1;

    unsigned int Nzp = Nz/2+1;
    size_t align = sizeof(Complex);

    // in
    array3<double> dens_r(Nx,Ny,Nz,align);
    array3<double> ergy_r(Nx,Ny,Nz,align);

    // out
    array3<Complex> dens_c(Nx,Ny,Nzp,align);
    array3<Complex> ergy_c(Nx,Ny,Nzp,align);

    // process
    rcfft3d dens_forward(Nx,Ny,Nz,dens_r,dens_c);
    rcfft3d ergy_forward(Nx,Ny,Nz,ergy_r,ergy_c);

    // ---------------------------------------------------------------------- //
    // Fill in real valued input arrays

    for (uint i = 0 ; i < Nx ; i++)
    {
        double x = (i + 0.5)*dx;

        for (uint j = 0 ; j < Ny ; j++)
        {
            double y = (j + 0.5)*dy;

            for (uint k = 0 ; k < Nz ; k++)
            {
                double z = (k + 0.5)*dz;

                dvec pos = {x,y,z};
                dvec bdata;

                uint bindex, cindex;
                meshinfo.get_cell_index(pos,bindex,cindex);

                ddens.get_block_data(bindex, bdata);
                double dens = bdata[cindex];

                dvelx.get_block_data(bindex, bdata);
                double velx = bdata[cindex];

                dvely.get_block_data(bindex, bdata);
                double vely = bdata[cindex];

                dvelz.get_block_data(bindex, bdata);
                double velz = bdata[cindex];

                dens_r(i,j,k) = dens;
                ergy_r(i,j,k) = norm3d(velx,vely,velz);

            }
        }
        cerr << i << endl;
    }

    // ---------------------------------------------------------------------- //
    // Compute Fourier Spectra

    dens_forward.fft0Normalized(dens_r,dens_c);
    ergy_forward.fft0Normalized(ergy_r,ergy_c);

    vector<double> dens_avg(Nx);
    vector<double> ergy_avg(Nx);
    vector<uint> avg_n(Nx);


    // ---------------------------------------------------------------------- //
    // Calculate spherical shell average

    for (uint i = 0 ; i < Nx ; i++)
    {
        for (uint j = 0 ; j < Ny ; j++)
        {
            for (uint k = 0 ; k < Nzp ; k++)
            {
                int p = i - Nx/2;
                int q = j - Ny/2;
                int r = k;

                double K = norm3d(p,q,r);
                dens_avg[round(K)] += abs(dens_c(i,j,k));
                ergy_avg[round(K)] += abs(ergy_c(i,j,k));
                avg_n[round(K)]++;
            }
        }
    }

    // ---------------------------------------------------------------------- //
    // Calculate power spectrum

    for (uint K = 0 ; K < Nx ; K++)
    {
        if ( avg_n[K] > 0 )
            cout
                << K << "\t"
                << K*K * dens_avg[K]/avg_n[K] << "\t"
                << K*K * ergy_avg[K]/avg_n[K] << "\t"
                << endl;
        else
            cout << K << 0 << 0 << endl;
    }

    return 0;

    // ---------------------------------------------------------------------- //
    // Spit out results

    for (uint i = 0 ; i < Nx ; i++)
    {
        double x = (i + 0.5)*dx;

        for (uint j = 0 ; j < Ny ; j++)
        {
            double y = (j + 0.5)*dy;

            for (uint k = 0 ; k < Nzp ; k++)
            {
                double z = (k + 0.5)*dz;

                cout
                    << i << "\t"
                    << j << "\t"
                    << k << "\t"

                    << x << "\t"
                    << y << "\t"
                    << z << "\t"

                    << abs( dens_c(i,j,k) ) << "\t"
                    << arg( dens_c(i,j,k) ) << "\t"

                    << abs( ergy_c(i,j,k) ) << "\t"
                    << arg( ergy_c(i,j,k) ) << "\t"

                    << endl;
            }
        }
        cout << endl;
    }

    return 0;
}
