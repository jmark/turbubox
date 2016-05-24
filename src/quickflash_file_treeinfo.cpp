// C++ program file quickflash_file_treeinfo.cpp

/*
  By Nathan C. Hearn
     October 4, 2008

  Reader for block centers, widths, and topology from Flash data files.
*/


#include "quickflash_file_treeinfo.hpp"
#include <hdf5.h>
#include "quickflash_except.hpp"
#include "quickflash_hdf5.hpp"
#include "quickflash_mesh_boundarydefs.hpp"
#include "quickflash_file_flashdefs.hpp"


namespace QuickFlash
{
namespace File
{

// Class TreeInfo

TreeInfo::TreeInfo() : 
  dims(0), num_blocks(0), block_center(), block_width(), block_parent_index(), 
  block_neighbor_indexes_lo(), block_neighbor_indexes_hi(), 
  block_neighbor_types_lo(), block_neighbor_types_hi(), 
  block_boundary_types_lo(), block_boundary_types_hi(), 
  block_child_indexes(), block_type(),
  block_refine_level(), min_refine_level(0), max_refine_level(0), 
  block_process_id()
  { reset(); }


TreeInfo::TreeInfo(const hid_t file_id, const FlashVersion flash_version, 
		   const unsigned int mesh_dims) : 
  dims(0), num_blocks(0), block_center(), block_width(), block_parent_index(), 
  block_neighbor_indexes_lo(), block_neighbor_indexes_hi(), 
  block_neighbor_types_lo(), block_neighbor_types_hi(), 
  block_boundary_types_lo(), block_boundary_types_hi(), 
  block_child_indexes(), block_type(),
  block_refine_level(), min_refine_level(0), max_refine_level(0),
  block_process_id()
  { reset(file_id, flash_version, mesh_dims); }


TreeInfo::TreeInfo(const TreeInfo & source) : 
  dims(0), num_blocks(0), block_center(), block_width(), block_parent_index(), 
  block_neighbor_indexes_lo(), block_neighbor_indexes_hi(), 
  block_neighbor_types_lo(), block_neighbor_types_hi(), 
  block_boundary_types_lo(), block_boundary_types_hi(), 
  block_child_indexes(), block_type(),
  block_refine_level(), min_refine_level(0), max_refine_level(0),
  block_process_id()
  { reset(source); }


void TreeInfo::reset()
  {
  dims = 0;

  num_blocks = 0;

  block_center.clear();
  block_width.clear();

  block_parent_index.clear();

  block_neighbor_indexes_lo.clear();
  block_neighbor_indexes_hi.clear();

  block_neighbor_types_lo.clear();
  block_neighbor_types_hi.clear();

  block_boundary_types_lo.clear();
  block_boundary_types_hi.clear();

  block_child_indexes.clear();

  block_type.clear();

  block_refine_level.clear();
  min_refine_level = 0;
  max_refine_level = 0;

  block_process_id.clear();
  }


void TreeInfo::reset(const hid_t file_id, const FlashVersion flash_version, 
		     const unsigned int mesh_dims)
  {
  reset();

  if (file_id >= 0)
    {
    read_block_info(file_id, mesh_dims);
    read_neighbor_info(flash_version, file_id);
    read_refinement_levels(file_id);
    read_block_types(file_id);
    read_process_ids(file_id);
    }
  else
    throw Except("Improper HDF5 file handle", __FILE__, __LINE__);
  }


void TreeInfo::reset(const TreeInfo & source)
  {
  if (&source != this)
    {
    dims = source.dims;

    num_blocks = source.num_blocks;

    block_center = source.block_center;
    block_width = source.block_width;

    block_parent_index = source.block_parent_index;

    block_neighbor_indexes_lo = source.block_neighbor_indexes_lo;
    block_neighbor_indexes_hi = source.block_neighbor_indexes_hi;

    block_neighbor_types_lo = source.block_neighbor_types_lo;
    block_neighbor_types_hi = source.block_neighbor_types_hi;

    block_boundary_types_lo = source.block_boundary_types_lo;
    block_boundary_types_hi = source.block_boundary_types_hi;

    block_child_indexes = source.block_child_indexes;

    block_type = source.block_type;

    block_refine_level = source.block_refine_level;
    min_refine_level = source.min_refine_level;
    max_refine_level = source.max_refine_level;

    block_process_id = source.block_process_id;
    }
  }


void TreeInfo::read_block_info(const hid_t file_id, 
			       const unsigned int mesh_dims)
  {
  using HDF5::read_data;

  dims = mesh_dims;

  // Get the block centers -- force vectors to have proper dimensions

  if (read_data(file_id, Coordinates_Name, mesh_dims, block_center) < 0)
    throw Except("Unable to read block centers", __FILE__, __LINE__);

  num_blocks = block_center.size();

  if (num_blocks < 1)
    throw Except("No blocks present in file", __FILE__, __LINE__);

  // Get the block widths -- force vectors to have proper dimensions

  if (read_data(file_id, BlockWidth_Name, mesh_dims, block_width) < 0)
    throw Except("Unable to read block widths", __FILE__, __LINE__);
  }


void TreeInfo::read_neighbor_info(const FlashVersion flash_version,
				  const hid_t file_id)
  {
  using HDF5::read_data;

  // WARNING: Call read_block_info before calling this function!

  // Get the neighbor indexes

  block_parent_index.resize(num_blocks);

  block_neighbor_indexes_lo.resize(num_blocks);
  block_neighbor_indexes_hi.resize(num_blocks);

  block_neighbor_types_lo.resize(num_blocks);
  block_neighbor_types_hi.resize(num_blocks);

  block_boundary_types_lo.resize(num_blocks);
  block_boundary_types_hi.resize(num_blocks);

  block_child_indexes.resize(num_blocks);

  std::vector< std::vector<int> > neighbor_index_data;

  if (read_data(file_id, Tree_Struct_Name, neighbor_index_data) < 0)
    throw Except("Unable to read block neighbor data", __FILE__, __LINE__);

  if (neighbor_index_data.size() != num_blocks)
    throw Except("Incorrect number of neighbor index data sets", __FILE__, 
		 __LINE__);

  // The neighbor index data is arranged as: <faces> parent <children>

  const unsigned int neighbors_start = 0;  // "faces"
  const unsigned int num_neighbors = 2 * dims;

  const unsigned int parent_start = neighbors_start + num_neighbors;

  const unsigned int children_start = parent_start + 1;
  const unsigned int num_children = 1 << dims;

  const unsigned int neighbor_index_size = num_neighbors + 1 + num_children;

  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
    const std::vector<int> & index_vector = neighbor_index_data[block_index];

    if (index_vector.size() != neighbor_index_size)
      throw Except("Improper neighbor index vector length", __FILE__, 
		   __LINE__);

    // Get the parent (don't forget Fortran offset-by-one)

    const int parent = index_vector[parent_start] - 1;

    if (parent >= 0)
      block_parent_index[block_index] = parent;
    else
      block_parent_index[block_index] = block_index;  // Roots point to selves

    // Get the neighbors

    std::vector<unsigned int> & neighbor_indexes_lo
      = block_neighbor_indexes_lo[block_index];

    std::vector<unsigned int> & neighbor_indexes_hi
      = block_neighbor_indexes_hi[block_index];

    std::vector<Mesh::NeighborType> & neighbor_types_lo
      = block_neighbor_types_lo[block_index];

    std::vector<Mesh::NeighborType> & neighbor_types_hi
      = block_neighbor_types_hi[block_index];

    std::vector<Mesh::BoundaryType> & boundary_types_lo
      = block_boundary_types_lo[block_index];

    std::vector<Mesh::BoundaryType> & boundary_types_hi
      = block_boundary_types_hi[block_index];

    neighbor_indexes_lo.resize(dims);
    neighbor_indexes_hi.resize(dims);

    neighbor_types_lo.resize(dims);
    neighbor_types_hi.resize(dims);

    boundary_types_lo.resize(dims);
    boundary_types_hi.resize(dims);

    for (unsigned int axis_index = 0; axis_index < dims; axis_index++)
      {
      const unsigned int elem_index = neighbors_start + (2 * axis_index);

      // Get the neighbor blocks indexes (don't forget Fortran offset-by-one)

      const int lo_neighbor = index_vector[elem_index];
      const int hi_neighbor = index_vector[elem_index + 1];

      Mesh::NeighborType lo_neighbor_type = Mesh::NoNeighbor;
      Mesh::NeighborType hi_neighbor_type = Mesh::NoNeighbor;

      Mesh::BoundaryType lo_boundary_type = Mesh::NoBoundary;
      Mesh::BoundaryType hi_boundary_type = Mesh::NoBoundary;

      // Determine the neighbor types

      switch (flash_version)
	{
	case Flash2 :

	  if (Mesh::get_neighbor_type_flash2(lo_neighbor, 
					     lo_neighbor_type,
					     lo_boundary_type) < 0)
	    throw Except("Invalid neighbor code", __FILE__, __LINE__);

	  if (Mesh::get_neighbor_type_flash2(hi_neighbor, 
					     hi_neighbor_type,
					     hi_boundary_type) < 0)
	    throw Except("Invalid neighbor code", __FILE__, __LINE__);

	  break;

	case Flash3 :

	  if (Mesh::get_neighbor_type_flash3(lo_neighbor, 
					     lo_neighbor_type,
					     lo_boundary_type) < 0)
	    throw Except("Invalid neighbor code", __FILE__, __LINE__);

	  if (Mesh::get_neighbor_type_flash3(hi_neighbor, 
					     hi_neighbor_type,
					     hi_boundary_type) < 0)
	    throw Except("Invalid neighbor code", __FILE__, __LINE__);

	  break;

	default :

	  throw Except("Unrecognized Flash code version", __FILE__, __LINE__);
	}

      // Assign indexes (non-valid neighbors just point to the current block)

      unsigned int lo_neighbor_id = block_index;  // Default value
      unsigned int hi_neighbor_id = block_index;  // Default value

      if (lo_neighbor_type == Mesh::ValidNeighbor)
	lo_neighbor_id = static_cast<unsigned int>(lo_neighbor - 1);

      if (hi_neighbor_type == Mesh::ValidNeighbor)
	hi_neighbor_id = static_cast<unsigned int>(hi_neighbor - 1);

      // Save the data

      neighbor_indexes_lo[axis_index] = lo_neighbor_id;
      neighbor_indexes_hi[axis_index] = hi_neighbor_id;

      neighbor_types_lo[axis_index] = lo_neighbor_type;
      neighbor_types_hi[axis_index] = hi_neighbor_type;

      boundary_types_lo[axis_index] = lo_boundary_type;
      boundary_types_hi[axis_index] = hi_boundary_type;
      }

    // Count the number of valid children (don't forget Fortran offset-by-one)

    unsigned int child_count = 0;

    for (unsigned int child_index = 0; child_index < num_children; 
	 child_index++)
      {
      const int child = index_vector[child_index + children_start] - 1;

      const unsigned int child_unsigned = static_cast<unsigned int>(child);

      if ((child >= 0) && (child_unsigned != block_index))
	child_count++;
      }

    // Copy the child data (don't forget Fortran offset-by-one)

    std::vector<unsigned int> & child_list = block_child_indexes[block_index];

    child_list.resize(child_count);

    unsigned int child_ptr = 0;

    for (unsigned int child_index = 0; child_index < num_children;
	 child_index++)
      {
      const int child = index_vector[child_index + children_start] - 1;

      const unsigned int child_unsigned = static_cast<unsigned int>(child);

      if ((child >= 0) && (child_unsigned != block_index))
	child_list[child_ptr++] = child_unsigned;
      }
    }    
  }


void TreeInfo::read_refinement_levels(const hid_t file_id)
  {
  using HDF5::read_data;

  // WARNING: Call read_block_info before calling this function!

  // Determine the refinement level for the blocks

  std::vector<int> block_refine_int;

  if (read_data(file_id, RefineLevel_Name, block_refine_int) < 0)
    throw Except("Unable to read block refinement levels", __FILE__, __LINE__);

  if (block_refine_int.size() != num_blocks)
    throw Except("Incompatible number of blocks in refinement levels",
		 __FILE__, __LINE__);

  min_refine_level = static_cast<unsigned int>(block_refine_int[0]);
  max_refine_level = min_refine_level;

  block_refine_level.resize(num_blocks);

  for (unsigned int i = 0; i < num_blocks; i++)
    {
    const unsigned int current_block_level 
      = static_cast<unsigned int>(block_refine_int[i]);

    block_refine_level[i] = current_block_level;

    if (current_block_level < min_refine_level)
      min_refine_level = current_block_level;
    else if (current_block_level > max_refine_level)
      max_refine_level = current_block_level;
    }
  }


void TreeInfo::read_block_types(const hid_t file_id)
  {
  using HDF5::read_data;

  // WARNING: Call read_block_info before calling this function!

  // Get the (integer) block types

  std::vector<int> block_type_ids;

  if (read_data(file_id, BlockType_Name, block_type_ids) < 0)
    throw Except("Unable to read block node type", __FILE__, __LINE__);

  if (block_type_ids.size() != num_blocks)
    throw Except("Incompatible number of blocks for node type",
		 __FILE__, __LINE__);

  block_type.resize(num_blocks);

  for (unsigned int i = 0; i < num_blocks; i++)
    block_type[i] = static_cast<unsigned int>(block_type_ids[i]);
  }


void TreeInfo::read_process_ids(const hid_t file_id)
  {
  using HDF5::member_present;
  using HDF5::read_data;

  // WARNING: Call read_block_info before calling this function!

  if (member_present(file_id, ProcessID_Name))
    {
    std::vector<int> temp_process_ids;

    if (read_data(file_id, ProcessID_Name, temp_process_ids) < 0)
      throw Except("Unable to read process IDs", __FILE__, __LINE__);

    if (temp_process_ids.size() != num_blocks)
      throw Except("Incompatible number of blocks for process IDs", __FILE__,
		   __LINE__);

    block_process_id.resize(num_blocks);

    bool ids_ok = true;

    for (unsigned int m = 0; m < num_blocks; m++)
      {
      const int current_id = temp_process_ids[m];

      if (current_id < 0)
	ids_ok = false;

      block_process_id[m] = static_cast<unsigned int>(current_id);
      }
    
    if (!ids_ok)
      throw Except("Process IDs out of range", __FILE__, __LINE__);
    }
  else if (num_blocks < 2)
    {
    block_process_id.resize(num_blocks);

    // We are probably in uniform grid mode -- give a dummy value

    if (num_blocks > 0)
      block_process_id[0] = 0;
    }
  else
    throw Except("Process IDs not present");
  }


}  // End namespace File
}  // End namespace QuickFlash
