// C++ program file quickflash_block_utils.cpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Block handling utilities.
*/


#include "quickflash_block_utils.hpp"
#include "quickflash_except.hpp"
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
			 std::vector<double> & child_maxbounds)
  {
  const unsigned int dims = current_minbounds.size();

  // ERROR CHECKING???

  child_minbounds = current_minbounds;
  child_maxbounds = current_maxbounds;

  for (unsigned int i = 0; i < dims; i++)
    {
    const double axis_center 
      = 0.5 * (current_minbounds[i] + current_maxbounds[i]);

    const unsigned int bit_mask = 1 << i;

    if ((local_child_index & bit_mask) != 0)
      child_minbounds[i] = axis_center;
    else
      child_maxbounds[i] = axis_center;
    }
  }


void invert_cell_index(const std::vector<unsigned int> & block_dims,
		       const unsigned int axis,
		       const std::vector<unsigned int> & cell_index,
		       std::vector<unsigned int> & inverted_cell_index)
  {
  // ERROR CHECKING???

  // Reflect the coordinate across the midpoint of the specified axis

  const unsigned int max_axis_cell_index = block_dims[axis] - 1;

  const unsigned int axis_cell_index = cell_index[axis];

  inverted_cell_index = cell_index;

  inverted_cell_index[axis] = max_axis_cell_index - axis_cell_index;
  }


void invert_cell_index_in_place(const std::vector<unsigned int> & block_dims,
				const unsigned int axis,
				std::vector<unsigned int> & cell_index)
  {
  // ERROR CHECKING???

  // Reflect the coordinate across the midpoint of the specified axis

  const unsigned int max_axis_cell_index = block_dims[axis] - 1;

  const unsigned int axis_cell_index = cell_index[axis];

  cell_index[axis] = max_axis_cell_index - axis_cell_index;
  }


unsigned int invert_cell_index(const std::vector<unsigned int> & block_dims,
			       const unsigned int axis,
			       const unsigned int cell_index)
  {
  // ERROR CHECKING???

  // Find the desired component

  unsigned int divisor = 1;

  for (unsigned int i = 0; i < axis; i++)
    divisor *= block_dims[i];

  const unsigned int num_axis_cells = block_dims[axis];

  const unsigned int axis_comp = (cell_index / divisor) % num_axis_cells;

  // Subtract the old and add the new

  const unsigned int new_axis_comp = (num_axis_cells - 1) - axis_comp;

  const unsigned int new_cell_index 
    = cell_index + (divisor * (new_axis_comp - axis_comp));

  return new_cell_index;
  }


double get_volume_overlap(const std::vector<double> & primary_min_coords,
			  const std::vector<double> & primary_max_coords,
			  const std::vector<double> & secondary_min_coords,
			  const std::vector<double> & secondary_max_coords)
  {
  double overlap_volume = 0.0;  // Default value

  const unsigned int dims = primary_min_coords.size();

  if ((primary_max_coords.size() != dims) 
      || (secondary_min_coords.size() != dims)
      || (secondary_max_coords.size() != dims))
    throw Except("Incompatible boundary vectors", __FILE__, __LINE__);

  // Get overlap

  bool overlapping = true;

  double volume = 1.0;

  for (unsigned int i = 0; i < dims; i++)
    {
    const double primary_min = primary_min_coords[i];
    const double primary_max = primary_max_coords[i];

    const double secondary_min = secondary_min_coords[i];
    const double secondary_max = secondary_max_coords[i];

    const double overlap_min = (primary_min > secondary_min) 
      ? primary_min : secondary_min;

    const double overlap_max = (primary_max < secondary_max) 
        ? primary_max : secondary_max;

    if (!(overlap_min < overlap_max))
      overlapping = false;

    // Get the volume and center position of the overlap

    volume *= overlap_max - overlap_min;
    }

  if (overlapping)
    overlap_volume = volume;

  return overlap_volume;
  }


}  // End namespace Utils
}  // End namespace Block
}  // End namespace QuickFlash
