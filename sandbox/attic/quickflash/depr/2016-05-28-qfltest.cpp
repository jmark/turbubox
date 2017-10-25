# include <string>
# include <vector>
# include <iostream>
# include <stdlib.h>
# include <cmath>

# include "quickflash_file_datafile.hpp"

using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char * argv[])
{
    // ---------------------------------------------------------------------- //
    // Read in Commandline args

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " filename" << " dbname" << endl;
        return 1;
    }

    int argPtr = 1;
    const std::string fname = argv[argPtr++];
    const std::string dname = argv[argPtr++];

    // ---------------------------------------------------------------------- //
    // Open the file
    const QuickFlash::File::DataFile dfile(fname);

    // ---------------------------------------------------------------------- //
    // Open datasets
    const QuickFlash::File::Dataset &ddens = dfile.get_dataset(dname);

    return 0;
}
