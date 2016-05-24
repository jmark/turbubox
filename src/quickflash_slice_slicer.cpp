// C++ program file quickflash_slice_slicer.cpp

/*
  By Nathan C. Hearn
     February 11, 2008

  Slicing routines.
*/


#include "quickflash_slice_slicer.hpp"
#include <vector>
#include "quickflash_file_siminfo.hpp"
#include "quickflash_file_dataset.hpp"
#include "quickflash_slice_slicegrid.hpp"
#include "quickflash_except.hpp"
#include "quickflash_block_blockdata.hpp"
#include "quickflash_block_blockinfo.hpp"


namespace QuickFlash
{
namespace Slice
{

void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, SliceGrid & grid, 
		  const double empty_value)
  {
  const unsigned int space_dims = meshinfo.get_dims();

  const std::vector<unsigned int> subsample_dims(space_dims, 1);

  slice_simple(meshinfo, dataset, subsample_dims, grid, empty_value);
  }


void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, 
		  const unsigned int start_pixel_index, 
		  const unsigned int num_proc_pixels,
		  SliceGrid & grid, const double empty_value)
  {
  const unsigned int space_dims = meshinfo.get_dims();
  
  const std::vector<unsigned int> subsample_dims(space_dims, 1);

  slice_simple(meshinfo, dataset, subsample_dims, start_pixel_index,
	       num_proc_pixels, grid, empty_value);
  }


void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, 
		  const std::vector<unsigned int> & subsample_dims,
		  SliceGrid & grid, const double empty_value)
  {
  const unsigned int num_pixels = grid.get_num_pixels();

  slice_simple(meshinfo, dataset, subsample_dims, 0, num_pixels, grid,
	       empty_value);
  }


void slice_simple(const File::MeshInfo & meshinfo,
		  const File::Dataset & dataset, 
		  const std::vector<unsigned int> & subsample_dims,
		  const unsigned int start_pixel_index,
		  const unsigned int num_proc_pixels,
		  SliceGrid & grid, const double empty_value)
  {
  // NOTE: Does not alter data in grid outside of specified pixel range

  const unsigned int space_dims = meshinfo.get_dims();

  if (subsample_dims.size() != space_dims)
    throw Except("Improper subsample_dims size", __FILE__, __LINE__);

  if (grid.get_space_dims() != space_dims)
    throw Except("Incompatible slice grid", __FILE__, __LINE__);

  unsigned int num_subsamples = 1;

  for (unsigned int n = 0; n < space_dims; n++)
    num_subsamples *= subsample_dims[n];

  if (num_subsamples < 1)
    throw Except("Non-positive value for subsample count", __FILE__, __LINE__);

  const unsigned int end_pixel_index 
    = (start_pixel_index + num_proc_pixels) - 1;

  const unsigned int num_total_pixels = grid.get_num_pixels();

  if (end_pixel_index >= num_total_pixels)
    throw Except("Pixel index out of range", __FILE__, __LINE__);

  // Fill the mesh

  std::vector< std::vector<double> > sample_coords;

  Block::BlockData<double> blockdata;

  bool block_loaded = false;

  const Block::BlockInfo * blockinfo_ptr = 0;

  for (unsigned int index = start_pixel_index; index <= end_pixel_index; 
       index++)
    {
    grid.get_pixel_subsamples(index, subsample_dims, sample_coords);

    // CHECK sample_coords SIZE AGAINST num_subsamples???

    if (sample_coords.size() != num_subsamples)
      throw Except("Error computing subsample positions", __FILE__, __LINE__);

    double accum_value = 0.0;

    unsigned int pixel_subsample_count = 0;

    for (unsigned int subsamp = 0; subsamp < num_subsamples; subsamp++)
      {
      const std::vector<double> & coords = sample_coords[subsamp];

      bool skip_sample = false;

      bool load_new_block = true;

      if (block_loaded)
	if (blockinfo_ptr->in_block(coords))
	  {
	  load_new_block = false;
	  skip_sample = false;
	  }

      if (load_new_block)
	{
	if (meshinfo.in_bounds(coords))
	  {
	  const unsigned int block_index = meshinfo.get_block_index(coords);

	  blockinfo_ptr = &(meshinfo.get_block_info(block_index));
	
	  dataset.get_block_data(block_index, blockdata);

	  block_loaded = true;
	  }
	else
	  skip_sample = true;
	}

      if (!skip_sample)
	{
	const unsigned int cell_index 
	  = blockinfo_ptr->get_nearest_cell(coords);

	const double cell_data = blockdata[cell_index];

	accum_value += cell_data;

	pixel_subsample_count++;
	}
      }

    // Store a value only if samples have been collected

    if (pixel_subsample_count > 0)
      grid[index] = accum_value / static_cast<double>(pixel_subsample_count);
    else
      grid[index] = empty_value;
    }
  }


}  // End namespace Slice
}  // End namespace QuickFlash
