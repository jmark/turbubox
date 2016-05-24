# include <string>
# include <vector>
# include <iostream>
# include <stdlib.h>
# include <cmath>

# include "H5File.hpp"
# include "H5DataSet.hpp"
# include "H5DataSpace.hpp"
# include "H5DataType.hpp"

# include "quickflash_file_datafile.hpp"
# include "quickflash_block_blockinfo.hpp"
# include "quickflash_file_meshinfo.hpp"

# include "Array.h"

# include <unistd.h>
# include <libgen.h>

using std::cout;
using std::cerr;
using std::endl;

typedef unsigned int uint;
typedef std::vector<double> dvec;

double get_data(
    const QuickFlash::File::Dataset *ds,
    const uint bindex,
    const uint cindex
) {
    dvec bdata;
    ds->get_block_data(bindex, bdata);
    return bdata[cindex];
}

template <typename T>
void write_ds (
    HighFive::File &file,
    const std::string &dname,
    const std::vector<T> &vec
) {
    file.createDataSet<T>(dname, HighFive::DataSpace::From(vec))
        .write(const_cast<std::vector<T>&>(vec));
}

template <typename T>
void write_ds (
    HighFive::File &file,
    const std::string &dname,
    const Array::array3<T> &arr
) {
    const HighFive::DataSpace Nxyz (arr.Nx()*arr.Ny()*arr.Nz());
    double *data = &*arr;
    file.createDataSet<T>(dname, Nxyz).write( data );
}

int main(int argc, char * argv[])
{
    // ---------------------------------------------------------------------- //
    // Read in Commandline args

    if (argc < 4)
    {
        cerr    << "Usage: " << argv[0] 
                << " <input file>" 
                << " <output file>" 
                << " <dbname> <dbname> ..." 
                << endl;
        return 1;
    }

    int argPtr = 1;
    const std::string qfla_file_name = argv[argPtr++];
    const std::string hdf5_file_name = argv[argPtr++];

    std::vector<std::string> dbnames;
    while (argPtr < argc)
    {
        dbnames.push_back(argv[argPtr++]);
    }

    const uint nrDbs = argc - 3;

    // use f*cking cache, otherwise it is f*cking slooow...
    const uint buffer_size = 20; // in number of blocks
    const uint cache_size = 20;  // in number of blocks

    // ---------------------------------------------------------------------- //
    // Open input file
    const QuickFlash::File::DataFile dfile(qfla_file_name, buffer_size );

    // Open output file
    HighFive::File h5file(hdf5_file_name,
          HighFive::File::ReadWrite 
        | HighFive::File::Create 
        | HighFive::File::Truncate );

    // ---------------------------------------------------------------------- //
    // Get mesh info
    const QuickFlash::File::MeshInfo &meshinfo = dfile.get_mesh_info();

    // ---------------------------------------------------------------------- //
    // Get physical box volume (uniform grid) in all directions
    const dvec &bmin = meshinfo.get_volume_minbounds();
    const dvec &bmax = meshinfo.get_volume_maxbounds();

    // ---------------------------------------------------------------------- //
    // Get smallest physical cell volume (uniform grid) in all directions

    const uint nblocks = meshinfo.get_num_blocks();
    uint maxrix = 0; // index of cell with maximal refine level
    uint rlevel = 0; // to be found maximal refine level

    // search for block with maximal refine level -> smallest cell volume
    for (uint bindex = 0 ; bindex < nblocks ; ++bindex)
    {
        const uint rlevel_ = meshinfo.get_refine_level(bindex);
        if (rlevel_ > rlevel)
        {
            rlevel = rlevel_;
            maxrix = bindex;
        }
    }

    const QuickFlash::Block::BlockInfo &binfo = meshinfo.get_block_info(maxrix);
    const dvec &cvol = binfo.get_cell_width();

    // ---------------------------------------------------------------------- //
    // Get number of cells (assuming uniform grid) in all directions

    // Box Volume
    const dvec &bvol = meshinfo.get_volume_width();

    // Vector with cell count in x,y and z direction
    const std::vector<uint> dims = {
            (uint)std::round(bvol[0]/cvol[0]),
            (uint)std::round(bvol[1]/cvol[1]),
            (uint)std::round(bvol[2]/cvol[2])
        };

    // ---------------------------------------------------------------------- //
    // Initialize arrays and databases

    std::vector < Array::array3<double>* > arrv;
    std::vector < const QuickFlash::File::Dataset* > dbsv;

    for ( std::string &dbname : dbnames )
    {
        auto *arr = new Array::array3<double>(dims[0],dims[1],dims[2]);
        auto &dbs = dfile.get_dataset(dbname);

        dbs.set_cache_size( cache_size );

        arrv.push_back(  arr );
        dbsv.push_back( &dbs );
    }

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

                uint bindex, cindex;
                meshinfo.get_cell_index(pos,bindex,cindex);

                for ( uint I = 0; I < nrDbs; ++I ) {
                    // DANGER: dbsv[I] and arrv[I] could be invalid pointers!
                    arrv[I]->operator()(i,j,k) = get_data(dbsv[I], bindex, cindex);
                }
            }
        }
    }

    // ---------------------------------------------------------------------- //
    // write to HDF5 file

    // META DATA
    write_ds(h5file, "DIMS"       ,dims);
    write_ds(h5file, "CELL_VOL"   ,cvol);
    write_ds(h5file, "BOX_VOL"    ,bvol);
    write_ds(h5file, "MIN_BOUNDS" ,bmin);
    write_ds(h5file, "MAX_BOUNDS" ,bmax);

    // PAYLOAD
    for ( uint i = 0;  i < dbnames.size(); ++i )
        write_ds(h5file, dbnames[i], *arrv[i] );

    return 0;
}
