// C++ header file quickflash_slice_slicer.hpp

/*
  By Nathan C. Hearn
     February 11, 2008

  Slicing routines.
*/


#ifndef QUICKFLASH_SLICE_SLICER_HPP
#define QUICKFLASH_SLICE_SLICER_HPP


#include <string>
#include <vector>
#include "quickflash_file_meshinfo.hpp"
#include "quickflash_file_dataset.hpp"
#include "quickflash_slice_slicegrid.hpp"


namespace QuickFlash
{
namespace Slice
{

void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, SliceGrid & grid,
		  const double empty_value=0.0);


void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, 
		  const unsigned int start_pixel_index, 
		  const unsigned int num_proc_pixels,
		  SliceGrid & grid, const double empty_value=0.0);


void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, 
		  const std::vector<unsigned int> & subsample_dims,
		  SliceGrid & grid, const double empty_value=0.0);


void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, 
		  const std::vector<unsigned int> & subsample_dims,
		  const unsigned int start_pixel_index,
		  const unsigned int num_proc_pixels,
		  SliceGrid & grid, const double empty_value=0.0);


}  // End namespace Slice
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SLICE_SLICER_HPP
