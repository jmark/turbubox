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

# include "ulz.hpp"

int main(int argc, char * argv[])
{
    // ---------------------------------------------------------------------- //
    // Read in Commandline args

    if (argc < 3)
    {
        std::cerr
                << "Usage: " << argv[0] 
                << " <input file>" 
                << " <dbname> <dbname> ..." 
                << std::endl;
        return 1;
    }

    int argPtr = 1;
    const std::string input_file_name = argv[argPtr++];

    std::vector<std::string> dbnames;
    while (argPtr < argc)
        dbnames.push_back(argv[argPtr++]);

    const uint nrDbs = dbnames.size();

    // ---------------------------------------------------------------------- //
    // Open input
    const HighFive::File infile(input_file_name, HighFive::File::ReadOnly);
    
    // ---------------------------------------------------------------------- //
    // META DATA
    std::vector<uint> dims(3);
    read_ds(infile,"DIMS",dims);
    
    // ---------------------------------------------------------------------- //
    // Initialize input arrays

    std::vector < Array::array3<double>* > arrv;
    
    for ( std::string &dbname : dbnames )
    {
        auto *arr  = new Array::array3<double>(dims[0],dims[1],dims[2]);
        read_ds(infile,dbname,*arr);
        arrv.push_back(  arr );
    }

    dvec AVG(nrDbs);
    const uint NNN = dims[0]*dims[1]*dims[2];

    // ---------------------------------------------------------------------- //
    // Calculate spherical shell cummulants

    nvec idx(3); // index vector
    dvec pos(3); // position vector

    for (idx[0] = 0 ; idx[0] < dims[0] ; ++idx[0])
    for (idx[1] = 0 ; idx[1] < dims[1] ; ++idx[1])
    for (idx[2] = 0 ; idx[2] < dims[2] ; ++idx[2])
    {
        for ( uint i = 0; i < nrDbs; ++i )
        {
            const double tmp = arrv[i]->operator()(idx[0],idx[1],idx[2]);
            AVG[i] += tmp/NNN;
        }
    }

    for ( uint I = 0; I < nrDbs; ++I )
        std::cout << AVG[I] << "\t";

    std::cout << std::endl;

    return 0;
}
