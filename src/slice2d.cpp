#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

# include "H5File.hpp"
# include "H5DataSet.hpp"
# include "H5DataSpace.hpp"
# include "H5DataType.hpp"

# include "Array.h"
# include "fftw++.h"

# include "ulz.hpp"

void usage()
{
    std::cerr
        << " Usage: " << "slice2d" 
        << " <input file>" 
        << " <dsname> <dsname> ..." 
        << " @" 
        << " <axis0> <axis1> <index>" 
        << std::endl;

    std::exit(1);
}

void die(const std::string msg)
{
    std::cerr << " " << msg << std::endl;
    std::exit(1); 
}

int main(int argc, char * argv[])
{
    // ---------------------------------------------------------------------- //
    // Read in Commandline args

    // least number of arguments
    if (argc < 7) usage();

    int argPtr = 1;
    const std::string input_file_name = argv[argPtr++];

    // delimiter used to distinguish the dataset name and the axis arguments
    const std::string delimiter = "@";

    // read in dataset names
    std::vector<std::string> dbnames;
    while (argPtr < argc && delimiter.compare(argv[argPtr]))
        dbnames.push_back(argv[argPtr++]);

    // It's bad when no dataset names were given.
    if(dbnames.empty()) usage();

    // dbnames size is needed later on.
    const uint nrDbs = dbnames.size();

    // We need exactly three more arguments.
    if(argc - ++argPtr != 3) usage();

    // vector with the axis information
    nvec axis(3);
    
    // store the two variying dimensions/axes
    axis[0] = std::stoul(argv[argPtr++]);
    axis[1] = std::stoul(argv[argPtr++]);

    if (axis[0] == axis[1]) die("axis0 and axis1 must not be the same!");
    if (axis[0] > 2 || axis[1] > 2) die("axis 0 and axis1 must be within 0 and 2!");

    // identify the constant dimension/axis
    axis[2] = 0 + 1 + 2 - axis[0] - axis[1];

    // store index for the constant dimension/axis
    uint component = std::stoul(argv[argPtr++]);

    // ---------------------------------------------------------------------- //
    // Open input file
    const HighFive::File infile(input_file_name, HighFive::File::ReadOnly);
    
    // ---------------------------------------------------------------------- //
    // META DATA
    nvec dims(3);
    dvec cvol(3);
    dvec bmin(3);

    read_ds(infile,"DIMS",dims);
    read_ds(infile,"CELL_VOL",cvol);
    read_ds(infile,"MIN_BOUNDS",bmin);

    if (component > dims[axis[2]]) die("Component not within bounds of array!");

    // ---------------------------------------------------------------------- //
    // Initialize input arrays

    std::vector < Array::array3<double>* > arrv;
    
    for ( std::string &dbname : dbnames )
    {
        auto *arr  = new Array::array3<double>(dims[0],dims[1],dims[2]);
        read_ds(infile,dbname,*arr);
        arrv.push_back(  arr );
    }

    // ---------------------------------------------------------------------- //

    nvec idx(3); // index vector
    dvec pos(3); // position vector

    //TODO: map axis to axis0, axis1 and (constant) axis2

    idx[axis[2]] = component;

    for (idx[axis[0]] = 0 ; idx[axis[0]] < dims[axis[0]] ; ++idx[axis[0]])
    {
        using std::cout;
        using std::endl;

        for (idx[axis[1]] = 0 ; idx[axis[1]] < dims[axis[1]] ; ++idx[axis[1]])
        {
            indexToPosition(bmin,cvol,idx,pos);                
           
            cout << "\t" << idx[axis[0]]; 
            cout << "\t" << idx[axis[1]]; 
                        
            cout << "\t" << pos[axis[0]]; 
            cout << "\t" << pos[axis[1]]; 

            for ( uint i = 0; i < nrDbs; ++i )
                cout << "\t" << arrv[i]->operator()(idx[0],idx[1],idx[2]);
            
            cout << endl;
        }
        cout << endl;
    }

    return 0;
}
