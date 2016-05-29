#pragma once

# include "H5File.hpp"
# include "Array.h"

typedef unsigned int uint;
typedef std::vector<uint> nvec;
typedef std::vector<double> dvec;

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
