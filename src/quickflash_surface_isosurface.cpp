// C++ program file quickflash_surface_isosurface.cpp

/*
  By Nathan C. Hearn
     July 16, 2008

  Class for encapsulating isosurface construction functions.
*/


#include "quickflash_surface_isosurface.hpp"
#include <list>
#include <vector>
#include "quickflash_file_meshinfo.hpp"
#include "quickflash_file_dataset.hpp"
#include "quickflash_surface_triangle.hpp"
#include "quickflash_surface_cube.hpp"
#include "quickflash_except.hpp"
#include "quickflash_surface_marchingtetrahedrons.hpp"


namespace QuickFlash
{
namespace Surface
{

// Class IsosurfaceGen

const unsigned int IsosurfaceGen::dims = 3;


void IsosurfaceGen::reset()
  {
  meshinfo = 0;
  dset = 0;

  limit_refine_level = false;
  max_refine_level = 0;

  limit_vol_min = false;
  vol_min_bounds.clear();

  limit_vol_max = false;
  vol_max_bounds.clear();
  }


void IsosurfaceGen::reset(const File::MeshInfo & mesh_info,
			  const File::Dataset & isosurface_dataset,
			  const bool limit_refinement_level,
			  const unsigned int max_refinement_level)
  {
  reset();

  if (mesh_info.get_dims() != dims)
    throw Except("Only three dimensional simulations supported",
			   __FILE__, __LINE__);

  meshinfo = &mesh_info;

  dset = &isosurface_dataset;

  set_refinement_limit(max_refinement_level, limit_refinement_level);
  }


void IsosurfaceGen::reset(const IsosurfaceGen & source)
  {
  if (&source != this)
    {
    meshinfo = source.meshinfo;
    dset = source.dset;

    limit_refine_level = source.limit_refine_level;
    max_refine_level = source.max_refine_level;

    limit_vol_min = source.limit_vol_min;
    vol_min_bounds = source.vol_min_bounds;

    limit_vol_max = source.limit_vol_max;
    vol_max_bounds = source.vol_max_bounds;
    }
  }


void IsosurfaceGen::set_refinement_limit(const unsigned int 
					   max_refinement_level,
					 const bool limit_refinement_level)
  {
  limit_refine_level = limit_refinement_level;

  if (limit_refine_level)
    {
    if (max_refinement_level < 1)
      throw Except("Flash refinement level must be at least one",
			     __FILE__, __LINE__);

    max_refine_level = max_refinement_level;
    }
  else
    max_refine_level = 0;
  }


void IsosurfaceGen::set_volume_limit(const bool limit_volume_min_bounds,
				     const std::vector<double> 
				       & volume_min_bounds,
				     const bool limit_volume_max_bounds,
				     const std::vector<double> 
				       & volume_max_bounds)
  {
  if (limit_volume_min_bounds)
    if (volume_min_bounds.size() != dims)
      throw Except("Volume min bounds must have three elements", __FILE__,
		   __LINE__);

  if (limit_volume_max_bounds)
    if (volume_max_bounds.size() != dims)
      throw Except("Volume max bounds must have three elements", __FILE__,
		   __LINE__);

  if (limit_volume_min_bounds && limit_volume_max_bounds)
    {
    bool elems_ok = true;  // Unless min >= max

    for (unsigned int i = 0; i < dims; i++)
      if (!(volume_min_bounds[i] < volume_max_bounds[i]))
	elems_ok = false;

    if (!elems_ok)
      throw Except("Volume min bounds must be less than max bounds", __FILE__,
		   __LINE__);
    }

  if (limit_volume_min_bounds)
    {
    limit_vol_min = true;
    vol_min_bounds = volume_min_bounds;
    }
  else
    {
    limit_vol_min = false;
    vol_min_bounds.clear();
    }

  if (limit_volume_max_bounds)
    {
    limit_vol_max = true;
    vol_max_bounds = volume_max_bounds;
    }
  else
    {
    limit_vol_max = false;
    vol_max_bounds.clear();
    }
  }


void IsosurfaceGen::add_block_surface(const unsigned int block_index,
				      const double isolevel_value,
				      std::list<Triangle> & triangle_list)
  {
  if ((meshinfo == 0) || (dset == 0))
    throw Except("Isosurface generator not initialized", __FILE__, __LINE__);

  const Block::BlockInfo & block_info = meshinfo->get_block_info(block_index);

  bool compute_surface = true;

  if (limit_vol_min || limit_vol_max)
    {
    const std::vector<double> & block_min_coords = block_info.get_min_coords();
    const std::vector<double> & block_max_coords = block_info.get_max_coords();

    if (!block_in_bounds(block_min_coords, block_max_coords))
      compute_surface = false;
    }

  if (compute_surface)
    {
    const std::vector<unsigned int> & block_dims = block_info.get_block_dims();

    const std::vector<double> & cell_width = block_info.get_cell_width();

    std::vector<double> buffer_width = cell_width;

    for (unsigned int i = 0; i < dims; i++)
      buffer_width[i] *= 0.56250;  // 9/16 of a cell width

    // Run through the data

    VectorCounter cube_counter(block_dims);

    const std::vector<unsigned int> corner_counter_dims(dims, 2);

    VectorCounter corner_counter(corner_counter_dims);

    std::vector<unsigned int> current_cell_index;
    std::vector<double> cell_pos;
    double cell_value = 0.0;

    Block::BlockData<double> block_data;

    dset->get_block_data(block_index, block_data);

    std::vector<unsigned int> temp_cell_index;
    std::vector<double> temp_cell_pos;

    bool neighbor_block_loaded = false;
    unsigned int prev_neighbor_block_index = 0;

    Block::BlockData<double> neighbor_block_data;

    CubeData<Triangle::DoubleVector> 
      vertex_positions;

    CubeData<double> vertex_data;

    while (cube_counter.in_bounds())
      {
      const std::vector<unsigned int> & start_cell_index 
	= cube_counter.get_state();

      corner_counter.reset_counter();

      bool skip_cube = false;

      while (corner_counter.in_bounds() && (!skip_cube))
	{
	const std::vector<unsigned int> & current_corner_index 
	  = corner_counter.get_state();

	current_cell_index = start_cell_index;

	for (unsigned int i = 0; i < dims; i++)
	  current_cell_index[i] += current_corner_index[i];

	bool in_main_block = true;

	for (unsigned int i = 0; i < dims; i++)
	  if (current_cell_index[i] >= block_dims[i])
	    in_main_block = false;

	if (in_main_block)
	  {
	  block_info.get_cell_center(current_cell_index, cell_pos);
	  cell_value = block_data[current_cell_index];
	  }
	else
	  {
	  // Get the nearest local cell

	  temp_cell_index = current_cell_index;

	  for (unsigned int i = 0; i < dims; i++)
	    {
	    const unsigned int block_max_axis_coord = block_dims[i] - 1;

	    if (temp_cell_index[i] > block_max_axis_coord)
	      temp_cell_index[i] = block_max_axis_coord;
	    }

	  block_info.get_cell_center(temp_cell_index, cell_pos);
	  
	  // Add (1/2 + 1/16) of the cell width along each axis that is outside

	  temp_cell_pos = cell_pos;

	  for (unsigned int i = 0; i < dims; i++)
	    if (current_cell_index[i] >= block_dims[i])
	      temp_cell_pos[i] += buffer_width[i];

	  // Check if this posiion is outside of the volume

	  if (!(meshinfo->in_bounds(temp_cell_pos)))
	    skip_cube = true;
	  else
	    {
	    // Get the block and cell index for this position

	    // NOTE: Max level for MeshInfo::get_block_index is one less than
	    //       the refinement level reported by Flash

	    const unsigned int neighbor_block_index 
	      = limit_refine_level 
	        ? meshinfo->get_block_index(temp_cell_pos, 
					    max_refine_level - 1)
	        : meshinfo->get_block_index(temp_cell_pos);

	    const Block::BlockInfo & neighbor_block_info 
	      = meshinfo->get_block_info(neighbor_block_index);

	    const unsigned int neighbor_cell_index 
	      = neighbor_block_info.get_nearest_cell(temp_cell_pos);

	    // Get the position and data for this cell

	    neighbor_block_info.get_cell_center(neighbor_cell_index, 
						temp_cell_pos);

	    bool load_new_block = true;

	    if (neighbor_block_loaded)
	      if (neighbor_block_index == prev_neighbor_block_index)
		load_new_block = false;

	    if (load_new_block)
	      {
	      dset->get_block_data(neighbor_block_index, neighbor_block_data);

	      prev_neighbor_block_index = neighbor_block_index;
	      neighbor_block_loaded = true;
	      }

	    cell_value = neighbor_block_data[neighbor_cell_index];

	    // Adjust the position of the vertex to line up with neighbor cell

	    for (unsigned int i = 0; i < dims; i++)
	      if (current_cell_index[i] >= block_dims[i])
		cell_pos[i] = temp_cell_pos[i];
	    }
	  }

	if (!skip_cube)
	  {
	  vertex_positions.get_elem(current_corner_index[0], 
				    current_corner_index[1], 
				    current_corner_index[2]) = cell_pos;

	  vertex_data.get_elem(current_corner_index[0], 
			       current_corner_index[1], 
			       current_corner_index[2]) = cell_value;
	  
	  corner_counter.increment();
	  }
	}

      if (!skip_cube)
	if (interesting_cube(isolevel_value, vertex_data))
	  {
	  std::list<Triangle> cube_triangles;

	  add_cube_triangles(isolevel_value, vertex_positions, vertex_data, 
			     cube_triangles);

	  if (limit_vol_min || limit_vol_max)
	    {
	    // Copy only in-bounds triangles

	    std::list<Triangle>::iterator tri_iter 
	      = cube_triangles.begin();

	    while (tri_iter != cube_triangles.end())
	      {
	      std::list<Triangle>::iterator next_tri_iter = tri_iter;
	      ++next_tri_iter;
		      
	      if (triangle_in_bounds(*tri_iter))
		triangle_list.splice(triangle_list.end(), cube_triangles, 
				     tri_iter);

	      tri_iter = next_tri_iter;
	      }
	    }
	  else
	    triangle_list.splice(triangle_list.end(), cube_triangles);
	  }

      cube_counter.increment();
      }
    }
  }


bool IsosurfaceGen::block_in_bounds(const std::vector<double> & min_coords,
				    const std::vector<double> & max_coords) 
  const
  {
  bool ret_val = true;  // Unless found otherwise, this is the default

  if (min_coords.size() != dims)
    throw Except("Min coords must be a three-vector", __FILE__, __LINE__);

  if (max_coords.size() != dims)
    throw Except("Max coords must be a three-vector", __FILE__, __LINE__);

  // ADDITIONAL ERROR CHECKING???  (COMPATIBLE MIN/MAX COORDS ELEMENTS???)

  if (limit_vol_min)
    {
    for (unsigned int i = 0; i < dims; i++)
      {
      const double axis_max_coord = max_coords[i];

      if (!(axis_max_coord > vol_min_bounds[i]))
	ret_val = false;
      }
    }

  if (limit_vol_max)
    {
    for (unsigned int i = 0; i < dims; i++)
      {
      const double axis_min_coord = min_coords[i];

      if (!(axis_min_coord < vol_max_bounds[i]))
	ret_val = false;
      }
    }

  return ret_val;
  }


bool IsosurfaceGen::triangle_in_bounds(const Triangle & triangle) const
  {
  bool ret_val = true;  // In bounds unless found otherwise

  const unsigned int tri_dims = 3;  // For the local template

  if (limit_vol_min || limit_vol_max)
    {
    for (unsigned int vertex_index = 0; vertex_index < 3; vertex_index++)
      {
      const Geometry::Vector<double, tri_dims> & vert_coords 
	= triangle[vertex_index];

      if (limit_vol_min)
	for (unsigned int i = 0; i < tri_dims; i++)
	  if (!(vert_coords[i] > vol_min_bounds[i]))
	    ret_val = false;  // Out of bounds

      if (limit_vol_max)
	for (unsigned int i = 0; i < tri_dims; i++)
	  if (!(vert_coords[i] < vol_max_bounds[i]))
	    ret_val = false;  // Out of bounds
      }
    }

  return ret_val;
  }


}  // End namespace Surface
}  // End namespace QuickFlash
