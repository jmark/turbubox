// C++ header file quickflash_block_utils.hpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Utilities for blocks.
*/


#ifndef QUICKFLASH_BLOCK_UTILS_HPP
#define QUICKFLASH_BLOCK_UTILS_HPP


#include <vector>


namespace QuickFlash
{
namespace Block
{
namespace Utils
{

void get_subblock_bounds(const unsigned int local_child_index,
			 const std::vector<double> & current_minbounds,
			 const std::vector<double> & current_maxbounds,
			 std::vector<double> & child_minbounds,
			 std::vector<double> & child_maxbounds);


void invert_cell_index(const std::vector<unsigned int> & block_dims,
		       const unsigned int axis,
		       const std::vector<unsigned int> & cell_index,
		       std::vector<unsigned int> & inverted_cell_index);


void invert_cell_index_in_place(const std::vector<unsigned int> & block_dims,
				const unsigned int axis,
				std::vector<unsigned int> & cell_index);


unsigned int invert_cell_index(const std::vector<unsigned int> & block_dims,
			       const unsigned int axis,
			       const unsigned int cell_index);


double get_volume_overlap(const std::vector<double> & primary_min_coords,
			  const std::vector<double> & primary_max_coords,
			  const std::vector<double> & secondary_min_coords,
			  const std::vector<double> & secondary_max_coords);


}  // End namespace Utils
}  // End namespace Block
}  // End namespace QuickFlash


#endif  // QUICKFLASH_BLOCK_UTILS_HPP
