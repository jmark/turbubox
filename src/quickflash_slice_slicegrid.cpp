// C++ program file quickflash_slice_slicegrid.cpp

/*
  By Nathan C. Hearn
     February 11, 2008

  Grid for composing 2D slices.

  NOTE: Currently, SliceGrid only works with cartesian geometry.
*/


#include "quickflash_slice_slicegrid.hpp"
#include <vector>
#include <math.h>
#include "quickflash_except.hpp"
#include "quickflash_geometry.hpp"
#include "quickflash_counters.hpp"


namespace QuickFlash
{
namespace Slice
{

// Class SliceGrid

void SliceGrid::reset()
  {
  space_dims = 0;

  grid_info.reset();
  grid_data.reset();

  grid_center.clear();

  slice_width.clear();

  orthog_min = 0.0;
  orthog_width = 0.0;

  roll = 0.0;
  polar = 0.0;
  azim = 0.0;

  trans_mat.clear();
  }


void SliceGrid::reset(const std::vector<double> & center_coords, 
		      const std::vector<double> & slice_grid_width, 
		      const std::vector<unsigned int> & pixel_count, 
		      const double azimuthal_angle, const double polar_angle, 
		      const double roll_angle, const double orthogonal_width)
  {
  space_dims = center_coords.size();

  if (space_dims < 2)
    throw Except("At least two spatial dimensions required", __FILE__, 
		 __LINE__);

  if (slice_grid_width.size() != 2)
    throw Except("Grid width must have two dimensions", __FILE__, __LINE__);

  for (unsigned int m = 0; m < 2; m++)
    if (!(slice_grid_width[m] > 0.0))
      throw Except("Grid widths must be positive", __FILE__, __LINE__);

  if (pixel_count.size() != 2)
    throw Except("Pixel count must have two dimensions", __FILE__, __LINE__);

  for (unsigned int m = 0; m < 2; m++)
    if (pixel_count[m] < 1)
      throw Except("Pixel counts must be positive", __FILE__, __LINE__);

  if (orthogonal_width < 0.0)
    throw Except("Orthogonal width must be non-negative", __FILE__, __LINE__);

  // Set up the grid -- coordinates relative to center_coords

  std::vector<double> grid_min_coords(2);
  std::vector<double> grid_max_coords(2);

  for (unsigned int m = 0; m < 2; m++)
    {
    const double half_width = 0.5 * slice_grid_width[m];

    grid_min_coords[m] = -half_width;
    grid_max_coords[m] = half_width;
    }

  grid_info.reset(grid_min_coords, grid_max_coords, pixel_count, 
		  Geometry::Cartesian);
  
  grid_data.reset(pixel_count);

  const unsigned int num_cells = grid_data.get_num_cells();

  // Clear the data

  for (unsigned int index = 0; index < num_cells; index++)
    grid_data[index] = 0.0;

  // Set up the orthogonal domain width

  orthog_min = -0.5 * orthogonal_width;  // Starting coordinate
  orthog_width = orthogonal_width;

  // Set up the transformation data : for column vectors, T=azim*polar*roll

  roll = roll_angle;
  polar = polar_angle;
  azim = azimuthal_angle;

  grid_center = center_coords;

  slice_width = slice_grid_width;

  set_transform_matrix();
  }


void SliceGrid::reset(const SliceGrid & source)
  {
  if (&source != this)
    {
    reset_metadata(source);

    grid_data = source.grid_data;
    }
  }


void SliceGrid::reset_metadata(const SliceGrid & source)
  {
  // Copies meta data, but not the data itself

  if (&source != this)
    {
    reset();

    space_dims = source.space_dims;

    grid_info = source.grid_info;

    grid_data.reset(source.grid_data.get_block_dims());  // Uninitialized data

    grid_center = source.grid_center;

    slice_width = source.slice_width;

    orthog_min = source.orthog_min;
    orthog_width = source.orthog_width;

    roll = source.roll;
    polar = source.polar;
    azim = source.azim;

    trans_mat = source.trans_mat;
    }
  }


void SliceGrid::get_pixel_center(const unsigned int pixel_index, 
				 std::vector<double> & space_coords) const
  {
  std::vector<double> grid_pos;

  grid_info.get_cell_center(pixel_index, grid_pos);

  if (grid_pos.size() != 2)
    throw Except("Unable to determine pixel center", __FILE__, __LINE__);

  // Pad start_vect with zeros

  grid_pos.resize(space_dims, 0.0);

  transform_pos(grid_pos, space_coords);
  }


void SliceGrid::get_pixel_center(const std::vector<unsigned int> & pixel_index,
				 std::vector<double> & coords) const
  {
  const unsigned int scalar_index = grid_info.get_scalar_index(pixel_index);
  get_pixel_center(scalar_index, coords);
  }


void SliceGrid::get_pixel_subsamples(const unsigned int pixel_index, 
				 const std::vector<unsigned int> & sample_dims,
		     std::vector< std::vector<double> > & sample_space_coords)
  const
  {
  if (sample_dims.size() != space_dims)
    throw Except("Incompatible sample dimensions vector", __FILE__, __LINE__);

  VectorCounter sample_counter(sample_dims);

  const unsigned int num_samples = sample_counter.get_num_states();

  if (num_samples < 1)
    throw Except("Sample dims must contain positive values", __FILE__, 
		 __LINE__);

  std::vector<double> pixel_min_coords;
  std::vector<double> pixel_max_coords;

  grid_info.get_cell_bounds(pixel_index, pixel_min_coords, pixel_max_coords);

  // Determine sample spacing and initial sample center

  double start_z = 0.0;
  double delta_z = 0.0;

  if (space_dims > 2)
    {
    const unsigned num_z_samples = sample_dims[2];

    if (num_z_samples > 1)
      {
      start_z = orthog_min;
      delta_z = orthog_width / static_cast<double>(num_z_samples);
      }
    }

  std::vector<double> coord_deltas(space_dims, 0.0);

  for (unsigned int m = 0; m < 2; m++)
    {
    const double axis_width = pixel_max_coords[m] - pixel_min_coords[m];

    const unsigned int axis_samples = sample_dims[m];

    if (axis_samples > 1)
      coord_deltas[m] = axis_width / static_cast<double>(axis_samples);
    else
      coord_deltas[m] = axis_width;
    }

  if (space_dims > 2)
    coord_deltas[2] = delta_z;

  std::vector<double> start_coords(space_dims);

  for (unsigned int m = 0; m < 2; m++)
    start_coords[m] = pixel_min_coords[m];

  if (space_dims > 2)
    start_coords[2] = start_z;

  // Align start coordinates with center of first cell on grid

  for (unsigned int i = 0; i < space_dims; i++)
    start_coords[i] += 0.5 * coord_deltas[i];

  // Collect the samples

  sample_space_coords.resize(num_samples);

  std::vector<double> grid_coords(space_dims);

  unsigned int sample_index = 0;

  sample_counter.reset_counter();

  while (sample_counter.in_bounds())
    {
    const std::vector<unsigned int> & counter_state 
      = sample_counter.get_state();

    for (unsigned int i = 0; i < space_dims; i++)
      {
      const double axis_delta 
	= static_cast<double>(counter_state[i]) * coord_deltas[i];

      grid_coords[i] = start_coords[i] + axis_delta;
      }

    transform_pos(grid_coords, sample_space_coords[sample_index]);

    sample_index++;

    sample_counter.increment();
    }
  }


void SliceGrid::get_pixel_subsamples(const std::vector<unsigned int> 
				     & pixel_index, 
			   const std::vector<unsigned int> & sample_dims,
		     std::vector< std::vector<double> > & sample_space_coords)
  const
  {
  const unsigned int scalar_index = grid_info.get_scalar_index(pixel_index);
  get_pixel_subsamples(scalar_index, sample_dims, sample_space_coords);
  }


void SliceGrid::get_symmetry_pair(const unsigned int symmetry_axis, 
				  const unsigned int pixel_index,
				  double & symm_sum, double & symm_diff) const
  {
  // ERROR CHECKING ON SYMMETRY_AXIS, PIXEL_INDEX???

  const unsigned int symm_index 
    = grid_info.invert_cell_index(symmetry_axis, pixel_index);

  const double pixel_value = grid_data[pixel_index];

  const double symm_value = grid_data[symm_index];

  symm_sum = 0.5 * (pixel_value + symm_value);

  symm_diff = 0.5 * (pixel_value - symm_value);
  }


void SliceGrid::get_symmetry_pair(const unsigned int symmetry_axis, 
				const std::vector<unsigned int> & pixel_index,
				  double & symm_sum, double & symm_diff) const
  {
  // ERROR CHECKING ON SYMMETRY_AXIS, PIXEL_INDEX???

  const unsigned int scalar_index = grid_info.get_scalar_index(pixel_index);

  get_symmetry_pair(symmetry_axis, scalar_index, symm_sum, symm_diff);
  }


void SliceGrid::set_symmetry_sum(const unsigned int symmetry_axis, 
				 const SliceGrid & source)
  {
  if (&source == this)
    throw Except("Can not use self to store symmetry", __FILE__, __LINE__);

  reset_metadata(source);

  if (symmetry_axis >= grid_data.get_space_dims())
    throw Except("Symmetry axis index out of range", __FILE__, __LINE__);

  const unsigned int num_cells = grid_data.get_num_cells();

  for (unsigned int index = 0; index < num_cells; index++)
    {
    const unsigned int symm_index 
      = grid_info.invert_cell_index(symmetry_axis, index);

    const double pixel_value = source.grid_data[index];

    const double symm_value = source.grid_data[symm_index];

    grid_data[index] = 0.5 * (pixel_value + symm_value);
    }
  }
  

void SliceGrid::set_symmetry_diff(const unsigned int symmetry_axis, 
				  const SliceGrid & source)
  {
  if (&source == this)
    throw Except("Can not use self to store symmetry");

  reset_metadata(source);

  if (symmetry_axis >= grid_data.get_space_dims())
    throw Except("Symmetry axis index out of range");

  const unsigned int num_cells = grid_data.get_num_cells();

  for (unsigned int index = 0; index < num_cells; index++)
    {
    const unsigned int symm_index 
      = grid_info.invert_cell_index(symmetry_axis, index);

    const double pixel_value = source.grid_data[index];

    const double symm_value = source.grid_data[symm_index];

    grid_data[index] = 0.5 * (pixel_value - symm_value);
    }
  }
  


void SliceGrid::set_transform_matrix()
  {
  trans_mat.resize(space_dims);

  for (unsigned int i = 0; i < space_dims; i++)
    trans_mat[i].resize(space_dims);

  const double ca = cos(azim);
  const double sa = sin(azim);

  switch(space_dims)
    {
    case 2 :
      {
      trans_mat[0][0] = ca;
      trans_mat[0][1] = -sa;

      trans_mat[1][0] = sa;
      trans_mat[1][1] = ca;
      }

      break;

    case 3 :
      {
      // trans_mat -> 
      //      z_rotate(azimuthal) * y_rotate(polar_angle) * z_rotate(roll)

      const double cp = cos(polar);
      const double sp = sin(polar);

      const double cr = cos(roll);
      const double sr = sin(roll);

      trans_mat[0][0] = (ca * cp * cr) - (sa * sr);
      trans_mat[0][1] = -(ca * cp * sr) - (cr * sa);
      trans_mat[0][2] = ca * sp;

      trans_mat[1][0] = (ca * sr) + (cp * cr * sa);
      trans_mat[1][1] = (ca * cr) - (cp * sa * sr);
      trans_mat[1][2] = sa * sp;

      trans_mat[2][0] = -cr * sp;
      trans_mat[2][1] = sp * sr;
      trans_mat[2][2] = cp;
      }

      break;

    default :

      throw Except("Incompatible dimensions for transform matrix", __FILE__, 
		   __LINE__);
    }
  }


void SliceGrid::transform_pos(const std::vector<double> & grid_coords,
			      std::vector<double> & space_coords) const
  {
  const unsigned int grid_coords_dims = grid_coords.size();

  if (grid_coords_dims != space_dims)
    throw Except("Incompatible dimensions for grid_coords", __FILE__,
		 __LINE__);

  space_coords = grid_center;

  // Add on the rotated coordinates

  for (unsigned int i = 0; i < space_dims; i++)
    {
    const std::vector<double> & trans_row = trans_mat[i];

    double rot_coord = 0.0;

    for (unsigned int j = 0; j < space_dims; j++)
      rot_coord += grid_coords[j] * trans_row[j];

    space_coords[i] += rot_coord;
    }
  }


}  // End namespace Slice
}  // End namespace QuickFlash
