# include <string>
# include <vector>
# include <iostream>
# include <stdlib.h>
# include <cmath>

# include "quickflash_file_datafile.hpp"
# include "quickflash_block_blockinfo.hpp"
# include "quickflash_file_meshinfo.hpp"

# include "Array.h"

# include "H5File.hpp"
# include "H5DataSet.hpp"
# include "H5DataSpace.hpp"
# include "H5DataType.hpp"

# include <unistd.h>
# include <libgen.h>

using std::cout;
using std::cerr;
using std::endl;

typedef unsigned int uint;
typedef std::vector<double> dvec;

double get_data_point (
    const QuickFlash::File::Dataset &ds,
    const dvec &pos,
    const QuickFlash::File::MeshInfo &meshinfo
) {
    uint bindex, cindex;
    meshinfo.get_cell_index(pos,bindex,cindex);
    dvec bdata;
    ds.get_block_data(bindex, bdata);
    return bdata[cindex];
}

bool is_readable(const std::string &path)
{
    if ( -1 == access(path.data(), R_OK)) {
        return false;
    }
    return true;
}

bool is_writable(const std::string &path)
{
    char dir[200];
    size_t length = path.copy(dir,200);
    dir[length] = '\0';

    dirname(dir);

    if ( -1 == access(dir, W_OK)) {
        return false;
    }
    return true;
}


int main(int argc, char * argv[])
{
    // ---------------------------------------------------------------------- //
    // Read in Commandline args

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " filename" << endl;
        return 1;
    }

    int argPtr = 1;
    const std::string qfla_file_name = argv[argPtr++];
    const std::string hdf5_file_name = argv[argPtr++];

    if ( !is_readable(qfla_file_name) ) {
        cerr << "[ERROR] '" << qfla_file_name << "': no read access!" << endl;
        return 1;
    }

    if ( !is_writable(hdf5_file_name) ) {
        cerr << "[ERROR] '" << hdf5_file_name << "': no write access!" << endl;
        return 1;
    }

    // ---------------------------------------------------------------------- //
    // Open the file
    const QuickFlash::File::DataFile dfile(qfla_file_name);

    // ---------------------------------------------------------------------- //
    // Get mesh info
    const QuickFlash::File::MeshInfo meshinfo = dfile.get_mesh_info();

    // ---------------------------------------------------------------------- //
    // Get physical box volume (uniform grid) in all directions
    dvec bmin = meshinfo.get_volume_minbounds();
    dvec bmax = meshinfo.get_volume_maxbounds();

    // ---------------------------------------------------------------------- //
    // Get physical cell volume (uniform grid) in all directions
    // TODO: find "smallest" cell in AMR grid
    const QuickFlash::Block::BlockInfo binfo = meshinfo.get_block_info(bmin);
    dvec cvol = binfo.get_cell_width();

    // ---------------------------------------------------------------------- //
    // Get number of cells (uniform grid) in all directions

    // Box Volume
    dvec bvol = meshinfo.get_volume_width();

    // Vector with grid number in x,y and z direction
    std::vector<uint> dims = {
            (uint)std::round(bvol[0]/cvol[0]),
            (uint)std::round(bvol[1]/cvol[1]),
            (uint)std::round(bvol[2]/cvol[2])
        };

    // ---------------------------------------------------------------------- //
    // Initialize arrays

    Array::array3<double> adens (dims[0],dims[1],dims[2]);
    Array::array3<double> avelx (dims[0],dims[1],dims[2]);
    Array::array3<double> avely (dims[0],dims[1],dims[2]);
    Array::array3<double> avelz (dims[0],dims[1],dims[2]);

    {
        // ---------------------------------------------------------------------- //
        // Open datasets
        const QuickFlash::File::Dataset ddens = dfile.get_dataset("dens");
        const QuickFlash::File::Dataset dvelx = dfile.get_dataset("velx");
        const QuickFlash::File::Dataset dvely = dfile.get_dataset("vely");
        const QuickFlash::File::Dataset dvelz = dfile.get_dataset("velz");

        // ---------------------------------------------------------------------- //
        // Fill up arrays

        dvec pos(3);

        for (uint i = 0 ; i < dims[0] ; i++)
        {
            pos[0] = (i + 0.5)*cvol[0];
            cerr << "x = " << i << endl;

            for (uint j = 0 ; j < dims[1] ; j++)
            {
                pos[1] = (j + 0.5)*cvol[1];

                for (uint k = 0 ; k < dims[2] ; k++)
                {
                    pos[2] = (k + 0.5)*cvol[2];

                    adens(i,j,k) = get_data_point(ddens,pos,meshinfo);
                    avelx(i,j,k) = get_data_point(dvelx,pos,meshinfo);
                    avely(i,j,k) = get_data_point(dvely,pos,meshinfo);
                    avelz(i,j,k) = get_data_point(dvelz,pos,meshinfo);
                }
            }
        }

    }

    // ---------------------------------------------------------------------- //
    // write to HDF5 file

    {
        HighFive::File file(qfla_file_name,
                            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

        // META DATA
        file.createDataSet<uint>  ("DIMS"       , HighFive::DataSpace::From(dims)).write( dims );
        file.createDataSet<double>("CELL_VOL"   , HighFive::DataSpace::From(cvol)).write( cvol );
        file.createDataSet<double>("BOX_VOL"    , HighFive::DataSpace::From(bvol)).write( bvol );
        file.createDataSet<double>("MIN_BOUNDS" , HighFive::DataSpace::From(bmin)).write( bmin );
        file.createDataSet<double>("MAX_BOUNDS" , HighFive::DataSpace::From(bmax)).write( bmax );

        // PAYLOAD
        const HighFive::DataSpace Nxyz (dims[0]*dims[1]*dims[2]);

        { double *data = &*adens; file.createDataSet<double>("dens", Nxyz).write( data ); }
        { double *data = &*avelx; file.createDataSet<double>("velx", Nxyz).write( data ); }
        { double *data = &*avely; file.createDataSet<double>("vely", Nxyz).write( data ); }
        { double *data = &*avelz; file.createDataSet<double>("velz", Nxyz).write( data ); }
    }

    return 0;
}
