#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "quickflash_file_datafile.hpp"
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_file_meshinfo.hpp"

using std::cout;
using std::cerr;
using std::endl;

typedef unsigned int uint;
typedef std::vector<double> dvec;

const double xmin = 0.0;
const double xmax = 1.0;
const double ymin = 0.0;
const double ymax = 1.0;
const double zmin = 0.0;
const double zmax = 1.0;

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
    if (argc < 2)
    {
        cerr    << "Usage: " << argv[0] << " filename"
	            << endl;
        return 1;
    }

    // Read in Commandline args
    int argPtr = 1;
    const std::string file_name     = argv[argPtr++];

    // Open the file, Get mesh info, get datasets
    const QuickFlash::File::DataFile    & dfile(file_name);
    const QuickFlash::File::MeshInfo    & meshinfo  = dfile.get_mesh_info();
    const QuickFlash::File::Dataset     & ddens     = dfile.get_dataset("dens");

    // const uint nblocks = meshinfo.get_num_blocks();
    // const uint dims    = meshinfo.get_dims();

    // const std::vector<uint> block_dims        = meshinfo.get_block_dims();
    // const std::vector<uint> base_block_dims   = meshinfo.get_base_block_dims();

    const dvec min = meshinfo.get_volume_minbounds();
    const dvec max = meshinfo.get_volume_maxbounds();
    const dvec vol = meshinfo.get_volume_width();

    const QuickFlash::Block::BlockInfo & binfo = meshinfo.get_block_info(min);
    const dvec cvl = binfo.get_cell_width();

    const double dx = cvl[0];
    const double dy = cvl[1];
    const double dz = cvl[2];

    const uint Nx = std::round(vol[0]/cvl[0]);
    const uint Ny = std::round(vol[1]/cvl[1]);
    const uint Nz = std::round(vol[2]/cvl[2]);

    // Mesh Info output
    // cerr << "Dims [ "       << dims         << " ]" << endl;
    // cerr << "Num blocks [ " << nblocks   << " ]" << endl;
    // cerr << "Num block cells [ " << meshinfo.get_num_block_cells() << " ]" << endl;

    // cerr << "Block dimensions [ ";
    // for (uint i = 0; i < dims; i++)
    //     cerr << block_dims[i] << " ";
    // //cerr << "]" << endl;

    // cerr << "Base block dimensions [ ";
    // for (uint i = 0; i < dims; i++)
    //     cerr << base_block_dims[i] << " ";
    // cerr << "]" << endl;

    //cerr << "Min  Bounds [ " << min[0] << " " << min[1] << " " << min[2] << " ]" << endl;
    //cerr << "Max  Bounds [ " << max[0] << " " << max[1] << " " << max[2] << " ]" << endl;
    //cerr << "Vol  Widths [ " << vol[0] << " " << vol[1] << " " << vol[2] << " ]" << endl;
    //cerr << "Cell Widths [ " << cvl[0] << " " << cvl[1] << " " << cvl[2] << " ]" << endl;

    //ddens[2] = dfile.get_dataset("dens");
    //ddens[3] = dfile.get_dataset("dens");

    // const QuickFlash::File::Dataset & dvelx = dfile.get_dataset("velx");
    // const QuickFlash::File::Dataset & dvely = dfile.get_dataset("vely");
    // const QuickFlash::File::Dataset & dvelz = dfile.get_dataset("velz");

    // Open datasets
    //cerr << "Datasets loaded" << endl;


    for (uint i = 0 ; i < Nx ; i++)
    {
        double x = (i + 0.5)*dx;

        for (uint j = 0 ; j < Ny ; j++)
        {
            double y = (j + 0.5)*dy;
            double dens = 0.0;

            for (uint k = 0 ; k < Nz ; k++)
            {
                double z = (k + 0.5)*dz;
                dvec pos = {x,y,z};

                uint bindex, cindex;
                meshinfo.get_cell_index(pos,bindex,cindex);

                dvec bdata;
                ddens.get_block_data(bindex, bdata);
                dens += bdata[cindex];
            }
            cout
                    << i << "\t"
                    << j << "\t"
                    << x << "\t"
                    << y << "\t"
                    << dens << "\t"
                    << std::log(dens) << "\t"
                    << endl;
        }
        cout << endl;
    }

    return 0;
}
