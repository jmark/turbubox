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

typedef unsigned int uint;
typedef std::vector<double> dvec;
typedef std::vector<uint> nvec;

double norm3d (const double x, const double y, const double z)
{
    return std::sqrt( x*x + y*y + z*z );
}

template <typename T>
double norm3d (const std::vector<T> & v)
{
    return std::sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

template <typename T>
void read_ds
(
    const HighFive::File &file,
    const std::string &dname,
    std::vector<T> &vec
)
{
    const auto &dataset = file.getDataSet(dname);
    dataset.read(vec);
}

template <typename T>
void read_ds (
    const HighFive::File &file,
    const std::string &dname,
    Array::array3<T> &arr
)
{
    const auto &dataset = file.getDataSet(dname);

    std::vector<T> tmp;
    dataset.read(tmp);

    arr = tmp.data();
}

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
    // size of result vectors
    const uint NK = std::ceil(norm3d(dims[0]/2,dims[1]/2,dims[2]))+1;

    // ---------------------------------------------------------------------- //
    // Initialize input and output arrays

    // input vector of arrays
    std::vector < Array::array3<double>* > arrv;
    
    // output vector
    dvec avgd (NK);
    nvec avgn (NK);

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

                double K = norm3d(p,q,r);

                for ( uint I = 0; I < nrDbs; ++I )
                    avgd[std::round(K)] += K*K * arrv[I]->operator()(i,j,k);

                avgn[std::round(K)]++;
            }
        }
    }

    // ---------------------------------------------------------------------- //
    // average cummulations and print out

    for (uint K = 0 ; K < NK ; K++)
    {
        using std::cout;
        using std::endl;

        double avg = avgd[K];
        const uint count = avgn[K];

        avg = count > 0 ? avg / count : 0;

        cout << K << "\t";
        cout << avg << "\t";
        cout << endl;
    }

    return 0;
}
