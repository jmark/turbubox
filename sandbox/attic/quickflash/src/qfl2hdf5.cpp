# include <string>
# include <vector>
# include <iostream>
# include <stdlib.h>
# include <cmath>

# include "ulz.hpp"

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

double get_block_data_point (
    const QuickFlash::File::Dataset *ds,
    const uint bindex,
    const uint cindex
) {
    dvec bdata;
    ds->get_block_data(bindex, bdata);
    return bdata[cindex];
}

template <typename T>
void pretty_print(const std::vector<T> &v)
{
    using std::cerr;
    using std::endl;

    for (uint i = 0; i < 3; i++)
        cerr << v[i] << " ";
    cerr << "]" << endl;
}

int main(int argc, char * argv[])
{
    bool DEBUG = false;
    if(std::getenv("DEBUG_QFL2HDF5"))
        DEBUG = true;

    // ---------------------------------------------------------------------- //
    // Read in Commandline args

    if (argc < 4)
    {
        std::cerr    
                << "Usage: " << argv[0] 
                << " <input file>" 
                << " <output file>" 
                << " <dbname> <dbname> ..." 
                << std::endl;
        return 1;
    }

    int argPtr = 1;
    const std::string qfla_file_name = argv[argPtr++];
    const std::string hdf5_file_name = argv[argPtr++];

    std::vector<std::string> dbnames;
    while (argPtr < argc)
        dbnames.push_back(argv[argPtr++]);

    const uint nrDbs = dbnames.size();

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
    // Debug output

    if (DEBUG)
    {
        using std::cerr;

        cerr << "Box Bimensions [ ";
        pretty_print(dims);

        cerr << "Box Min Bounds [ ";
        pretty_print(bmin);

        cerr << "Box Max Bounds [ ";
        pretty_print(bmax);

        cerr << "Box Volume     [ ";
        pretty_print(bvol);

        cerr << "Cell Volume    [ ";
        pretty_print(cvol);
    }


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

    nvec idx(3); // index vector
    dvec pos(3); // position vector

    if (DEBUG)
        std::cerr << std::endl << "progress in i/x-direction: " << std::endl;

    for (idx[0] = 0 ; idx[0] < dims[0] ; ++idx[0])
    {
        if (DEBUG)
            std::cerr << idx[0] << std::endl;
        for (idx[1] = 0 ; idx[1] < dims[1] ; ++idx[1])
        {
            for (idx[2] = 0 ; idx[2] < dims[2] ; ++idx[2])
            {
                indexToPosition(bmin,cvol,idx,pos);                
                uint bindex, cindex;
                meshinfo.get_cell_index(pos,bindex,cindex);

                for ( uint i = 0; i < nrDbs; ++i )
                {
                    arrv[i]->operator()(idx[0],idx[1],idx[2]) 
                        = get_block_data_point(dbsv[i], bindex, cindex);
                }
            }
        }
    }

    // ---------------------------------------------------------------------- //
    // Get some runtime metadata like time, time step, etc.
    const QuickFlash::File::SimInfo &siminfo = dfile.get_sim_info();
    std::vector<double> siminfov(3);

    siminfov[0] = siminfo.get_sim_time();
    siminfov[1] = siminfo.get_sim_dt();
    siminfov[2] = siminfo.get_step_count();

    // ---------------------------------------------------------------------- //
    // write to HDF5 file

    // META DATA
    write_ds(h5file, "DIMS"       ,dims);
    write_ds(h5file, "CELL_VOL"   ,cvol);
    write_ds(h5file, "BOX_VOL"    ,bvol);
    write_ds(h5file, "MIN_BOUNDS" ,bmin);
    write_ds(h5file, "MAX_BOUNDS" ,bmax);
    write_ds(h5file, "SIM_INFO"   ,siminfov);

    // PAYLOAD
    for ( uint i = 0;  i < dbnames.size(); ++i )
        write_ds(h5file, dbnames[i], *arrv[i] );

    return 0;
}
