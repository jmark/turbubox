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

    if (argc < 2)
    {
        std::cerr
                << "Usage: " << argv[0] 
                << " <input file>" 
                << std::endl;
        return 1;
    }

    int argPtr = 1;
    const std::string input_file_name = argv[argPtr++];

    // ---------------------------------------------------------------------- //
    // Open input
    const HighFive::File infile(input_file_name, HighFive::File::ReadOnly);
    
    // ---------------------------------------------------------------------- //
    // META DATA
    std::vector<uint> dims(3);
    read_ds(infile,"DIMS",dims);
    
    // ---------------------------------------------------------------------- //
    // Initialize input arrays

    auto velx = Array::array3<double>(dims[0],dims[1],dims[2]);
    auto vely = Array::array3<double>(dims[0],dims[1],dims[2]);
    auto velz = Array::array3<double>(dims[0],dims[1],dims[2]);

    auto magx = Array::array3<double>(dims[0],dims[1],dims[2]);
    auto magy = Array::array3<double>(dims[0],dims[1],dims[2]);
    auto magz = Array::array3<double>(dims[0],dims[1],dims[2]);

    read_ds(infile,"velx",velx);
    read_ds(infile,"vely",vely);
    read_ds(infile,"velz",velz);

    read_ds(infile,"magx",magx);
    read_ds(infile,"magy",magy);
    read_ds(infile,"magz",magz);

    const uint NNN = dims[0]*dims[1]*dims[2];

    double avg_ekin(0);
    double avg_emag(0);

    // ---------------------------------------------------------------------- //

    nvec idx(3); // index vector
    dvec pos(3); // position vector

    for (idx[0] = 0 ; idx[0] < dims[0] ; ++idx[0])
    for (idx[1] = 0 ; idx[1] < dims[1] ; ++idx[1])
    for (idx[2] = 0 ; idx[2] < dims[2] ; ++idx[2])
    {
        const dvec ekin_tmp {
            velx(idx[0],idx[1],idx[2]),
            vely(idx[0],idx[1],idx[2]),
            velz(idx[0],idx[1],idx[2])
        };

        const dvec emag_tmp {
            magx(idx[0],idx[1],idx[2]),
            magy(idx[0],idx[1],idx[2]),
            magz(idx[0],idx[1],idx[2])
        };

        avg_ekin += norm3d(ekin_tmp)/NNN;
        avg_emag += norm3d(emag_tmp)/NNN;
    }

    std::cout << avg_ekin << "\t";
    std::cout << avg_emag << "\t";

    std::cout << std::endl;

    return 0;
}
