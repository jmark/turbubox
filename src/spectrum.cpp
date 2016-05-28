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

template <typename T>
void write_ds (
    HighFive::File &file,
    const std::string &dname,
    const Array::array3<T> &arr
) {
    const HighFive::DataSpace Nxyz (arr.Nx()*arr.Ny()*arr.Nz());
    T *data = &*arr;
    file.createDataSet<T>(dname, Nxyz).write( data );
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

int main(int argc, char * argv[])
{
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
    const std::string input_file_name = argv[argPtr++];
    const std::string output_file_name = argv[argPtr++];

    std::vector<std::string> dbnames;
    while (argPtr < argc)
        dbnames.push_back(argv[argPtr++]);

    // ---------------------------------------------------------------------- //
    // Open input and output files
    const HighFive::File infile(input_file_name, HighFive::File::ReadOnly);
    
    HighFive::File outfile(output_file_name, 
          HighFive::File::ReadWrite 
        | HighFive::File::Create 
        | HighFive::File::Truncate );

    // ---------------------------------------------------------------------- //
    // META DATA
    std::vector<uint> dims(3);
    read_ds(infile,"DIMS",dims);

    std::vector<uint> fdims(3);
    fdims[0] = dims[0];
    fdims[1] = dims[1];
    fdims[2] = dims[2]/2+1;

    write_ds(outfile, "DIMS" ,fdims);

    // ---------------------------------------------------------------------- //
    // compute fourier trafos and write to output file

    // fftw::maxthreads = get_max_threads();
    fftwpp::fftw::maxthreads = 1;

    for ( std::string &dbname : dbnames )
    {
        Array::array3<double>  real        (dims[0],dims[1],dims[2]);
        Array::array3<Complex> fourier     (fdims[0],fdims[1],fdims[2]);
        Array::array3<double>  fourier_abs (fdims[0],fdims[1],fdims[2]);
        Array::array3<double>  fourier_arg (fdims[0],fdims[1],fdims[2]);

        read_ds(infile,dbname,real);

        // prepare process
        fftwpp::rcfft3d fforward(dims[0],dims[1],dims[2],real,fourier);

        // compute fourier spectra
        fforward.fft0Normalized(real,fourier);

        // compute abs and arg of fourier
        for (uint i = 0 ; i < fdims[0] ; i++)
        for (uint j = 0 ; j < fdims[1] ; j++)
        for (uint k = 0 ; k < fdims[2] ; k++)
        {
            fourier_abs(i,j,k) = std::abs(fourier(i,j,k));
            fourier_arg(i,j,k) = std::arg(fourier(i,j,k));
        }

        write_ds(outfile, dbname + "_abs", fourier_abs);
        write_ds(outfile, dbname + "_arg", fourier_arg);
    }

    return 0;
}
