#pragma once

# include "H5File.hpp"
# include "Array.h"

typedef unsigned int uint;
typedef std::vector<int> zvec;
typedef std::vector<uint> nvec;
typedef std::vector<double> dvec;

// I would like to use this struct, but we stick to a
// 2D vector<double> for now.
//
// typedef struct {
//     double time;
//     double time_step;
//     // double step_count;
//     
// } SimInfo;

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
){
    const auto &dataset = file.getDataSet(dname);
    dataset.read(vec);
}

template <typename T>
void read_ds (
    const HighFive::File &file,
    const std::string &dname,
    Array::array3<T> &arr
){
    const auto &dataset = file.getDataSet(dname);

    std::vector<T> tmp;
    dataset.read(tmp);

    arr = tmp.data();
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
    T *data = &*arr;
    file.createDataSet<T>(dname, Nxyz).write( data );
}

void indexToPosition (
    const dvec &bmin, 
    const dvec &cvol, 
    const nvec &idx, 
    dvec &pos
){
    pos[0] = bmin[0] + cvol[0] * (idx[0] + 0.5);
    pos[1] = bmin[1] + cvol[1] * (idx[1] + 0.5);
    pos[2] = bmin[2] + cvol[2] * (idx[2] + 0.5);
}
