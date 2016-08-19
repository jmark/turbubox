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
        " Usage: slice"
        " <input file>" 
        " <dsname> <dsname> ..." 
        " @" 
        " (x|[0-9]+) (y|[0-9]+) (z|[0-9]+)" 
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

	// bounds of indecs
    nvec lower(3);
    nvec upper(3);
    std::vector<bool> axis(3,false); // false -> const, true -> variable
    std::vector<char> slice_args = {'x','y','z'};

	for (uint i = 0 ; i < 3 ; ++i)
	{
		const char *arg = argv[argPtr++];
		if (*arg == slice_args[i])
		{
			lower[i] = 0;
			upper[i] = dims[i];
            axis[i]  = true;
		}
		else
		{
			try
			{
				const uint idx = std::stoul(arg);
				lower[i] = idx;
				upper[i] = idx+1;

				if (upper[i] > dims[i]) 
					die("Given index not within bounds of the box!");
			}
			catch (const std::exception& e) {
				die("Could not parse a proper index number!");
			}
		}
	}
	
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

	using std::cout;
	using std::endl;

    for (idx[0] = lower[0] ; idx[0] < upper[0] ; ++idx[0])
	{
    for (idx[1] = lower[1] ; idx[1] < upper[1] ; ++idx[1])
    {
    for (idx[2] = lower[2] ; idx[2] < upper[2] ; ++idx[2])
    {
        indexToPosition(bmin,cvol,idx,pos);                
       
        cout << idx[0];
        cout << "\t"; 
        cout << idx[1]; 
        cout << "\t"; 
        cout << idx[2]; 
        cout << "\t"; 
                      
        cout << pos[0]; 
        cout << "\t"; 
        cout << pos[1]; 
        cout << "\t"; 
        cout << pos[2]; 
        cout << "\t"; 

        for ( uint k = 0; k < nrDbs; ++k )
        {
            cout << arrv[k]->operator()(idx[0],idx[1],idx[2]);
            cout << "\t"; 
        }
        
    if (axis[2]) cout << endl;
    }
    if (axis[1]) cout << endl;
    }
	if (axis[0]) cout << endl;
	}

    return 0;
}
