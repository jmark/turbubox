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
    // Open input and output files
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

    // ---------------------------------------------------------------------- //
    // Calculate spherical shell cummulants

    for (uint i = 0 ; i < dims[0] ; i++)
    {
        for (uint j = 0 ; j < dims[1] ; j++)
        {
            for (uint k = 0 ; k < dims[2] ; k++)
            {
                int p = i - dims[0]/2;
                int q = j - dims[1]/2;
                int r = k;

                double K = abs3d(p,q,r);

                std::cout << K << "\t";

                for ( uint I = 0; I < nrDbs; ++I )
                    std::cout << arrv[I]->operator()(i,j,k) << "\t";
                
                std::cout << std::endl;
            }
        }
    }

    //TODO: average it out, too

    return 0;
}
