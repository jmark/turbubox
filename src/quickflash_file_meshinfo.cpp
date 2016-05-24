// C++ program file quickflash_file_meshinfo.cpp

/*
  By Nathan C. Hearn
     October 4, 2008
*/


#include "quickflash_file_meshinfo.hpp"
#include <hdf5.h>
#include <vector>
#include "quickflash_except.hpp"
#include "quickflash_file_siminfo.hpp"
#include "quickflash_block_blocktype.hpp"
#include "quickflash_file_treeinfo.hpp"
#include "quickflash_tree_builder.hpp"
#include "quickflash_tree_block.hpp"
#include "quickflash_tree_nodedata.hpp"
#include "quickflash_geometry.hpp"
#include "quickflash_file_flashdefs.hpp"


namespace QuickFlash
{
namespace File
{

// Class MeshInfo

MeshInfo::MeshInfo() :
  mesh_info_set_up(false), valid_mesh(false), sim_time(0.0), dims(0), 
  num_blocks(0), num_leaf_blocks(0), 
  min_leaf_refine_level(0), max_leaf_refine_level(0),
  block_refine_level(), block_process_id(),
  base_block_dims(), 
  volume_minbounds(), volume_maxbounds(), volume_center(), volume_width(),
  block_dims(), num_block_cells(0), 
  tree_block(), node_ptrs(), block_type(), geom_type(Geometry::Cartesian),
  parent_indexes(), 
  neighbor_indexes_lo(), neighbor_indexes_hi(),
  neighbor_types_lo(), neighbor_types_hi(), 
  boundary_types_lo(), boundary_types_hi()
  { }


MeshInfo::MeshInfo(const hid_t file_id, const SimInfo & sim_info) :
  mesh_info_set_up(false), valid_mesh(false), sim_time(0.0), dims(0), 
  num_blocks(0), num_leaf_blocks(0), 
  min_leaf_refine_level(0), max_leaf_refine_level(0),
  block_refine_level(), block_process_id(),
  base_block_dims(), 
  volume_minbounds(), volume_maxbounds(), volume_center(), volume_width(),
  block_dims(), num_block_cells(0), 
  tree_block(), node_ptrs(), block_type(), geom_type(Geometry::Cartesian),
  parent_indexes(), 
  neighbor_indexes_lo(), neighbor_indexes_hi(),
  neighbor_types_lo(), neighbor_types_hi(), 
  boundary_types_lo(), boundary_types_hi()
  { reset(file_id, sim_info); }


MeshInfo::MeshInfo(const MeshInfo & source) :
  mesh_info_set_up(false), valid_mesh(false), sim_time(0.0), dims(0), 
  num_blocks(0), num_leaf_blocks(0), 
  min_leaf_refine_level(0), max_leaf_refine_level(0),
  block_refine_level(), block_process_id(),
  base_block_dims(), 
  volume_minbounds(), volume_maxbounds(), volume_center(), volume_width(),
  block_dims(), num_block_cells(0), 
  tree_block(), node_ptrs(), block_type(), geom_type(Geometry::Cartesian),
  parent_indexes(), 
  neighbor_indexes_lo(), neighbor_indexes_hi(),
  neighbor_types_lo(), neighbor_types_hi(), 
  boundary_types_lo(), boundary_types_hi()
  { reset(source); }


void MeshInfo::reset()
  {
  mesh_info_set_up = false;

  valid_mesh = false;

  sim_time = 0.0;

  dims = 0;

  num_blocks = 0;
  num_leaf_blocks = 0;

  min_leaf_refine_level = 0;
  max_leaf_refine_level = 0;

  block_refine_level.clear();
  block_process_id.clear();

  base_block_dims.clear();

  volume_minbounds.clear();
  volume_maxbounds.clear();

  volume_center.clear();
  volume_width.clear();

  block_dims.clear();
  num_block_cells = 0;

  tree_block.reset();

  node_ptrs.clear();

  block_type.clear();

  geom_type = Geometry::Cartesian;

  parent_indexes.clear();

  neighbor_indexes_lo.clear();
  neighbor_indexes_hi.clear();

  neighbor_types_lo.clear();
  neighbor_types_hi.clear();

  boundary_types_lo.clear();
  boundary_types_hi.clear();
  }

  
void MeshInfo::reset(const hid_t file_id, const SimInfo & sim_info)
  {
  reset();

  if (file_id < 0)
    throw Except("Improper HDF5 file ID", __FILE__, __LINE__);

  get_sim_params(sim_info);

  if (valid_mesh)
    construct_tree(file_id, sim_info);

  mesh_info_set_up = true;
  }


void MeshInfo::reset(const MeshInfo & source)
  {
  if (&source != this)
    {
    mesh_info_set_up = source.mesh_info_set_up;

    valid_mesh = source.valid_mesh;

    sim_time = source.sim_time;

    dims = source.dims;

    num_blocks = source.num_blocks;
    num_leaf_blocks = source.num_leaf_blocks;

    min_leaf_refine_level = source.min_leaf_refine_level;
    max_leaf_refine_level = source.max_leaf_refine_level;

    block_refine_level = source.block_refine_level;
    block_process_id = source.block_process_id;

    base_block_dims = source.base_block_dims;

    volume_minbounds = source.volume_minbounds;
    volume_maxbounds = source.volume_maxbounds;

    volume_center = source.volume_center;
    volume_width = source.volume_width;

    block_dims = source.block_dims;
    num_block_cells = source.num_block_cells;

    Tree::reset_tree_block(source.tree_block, tree_block);

    set_node_data_pointers();

    block_type = source.block_type;

    geom_type = source.geom_type;

    parent_indexes = source.parent_indexes;

    neighbor_indexes_lo = source.neighbor_indexes_lo;
    neighbor_indexes_hi = source.neighbor_indexes_hi;

    neighbor_types_lo = source.neighbor_types_lo;
    neighbor_types_hi = source.neighbor_types_hi;

    boundary_types_lo = source.boundary_types_lo;
    boundary_types_hi = source.boundary_types_hi;
    }
  }


bool MeshInfo::in_bounds(const std::vector<double> & position) const
  {
  bool ret_val = true;

  // ERROR CHECKING???

  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  for (unsigned int i = 0; i < dims; i++)
    {
    const double coord = position[i];

    if ((coord < volume_minbounds[i]) || (!(coord < volume_maxbounds[i])))
      ret_val = false;
    }

  return ret_val;
  }


void MeshInfo::get_tree_iter(const std::vector<double> & position,
			    Tree::ConstTreeIter & tree_iter) const
  {
  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  // Identify the tree for this position

  const unsigned int num_trees = tree_block.get_num_cells();

  bool found_tree = false;

  unsigned int tree_index = 0;

  while ((!found_tree) && (tree_index < num_trees))
    {
    if (tree_block[tree_index].get_data().get_block_info().in_block(position))
      found_tree = true;
    else
      tree_index++;
    }

  if (!found_tree)
    throw Except("Tree not found", __FILE__, __LINE__);

  Tree::ConstTreeNodeIter node_iter = tree_block[tree_index].get_node_iter();

  bool found_it = false;

  while ((!found_it) && (node_iter.valid_node()))
    {
    const Block::BlockInfo & node_info = node_iter.get_data().get_block_info();

    if (!(node_info.in_block(position)))
      node_iter.next_skipBranch();
    else
      {
      if (node_iter.children_present())
	node_iter.next();
      else
	found_it = true;
      }
    }

  if (!found_it)
    throw Except("Block not found", __FILE__, __LINE__);

  const unsigned int block_index = node_iter.get_data().get_block_index();

  tree_iter = node_ptrs[block_index];
  }


void MeshInfo::get_tree_iter(const std::vector<double> & position,
			    const unsigned int max_level,
			    Tree::ConstTreeIter & tree_iter) const
  {
  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  // Identify the tree for this position

  const unsigned int num_trees = tree_block.get_num_cells();

  bool found_tree = false;

  unsigned int tree_index = 0;

  while ((!found_tree) && (tree_index < num_trees))
    {
    if (tree_block[tree_index].get_data().get_block_info().in_block(position))
      found_tree = true;
    else
      tree_index++;
    }

  if (!found_tree)
    throw Except("Tree not found", __FILE__, __LINE__);

  Tree::ConstTreeNodeIter node_iter = tree_block[tree_index].get_node_iter();

  bool found_it = false;

  while ((!found_it) && (node_iter.valid_node()))
    {
    const Block::BlockInfo & node_info = node_iter.get_data().get_block_info();

    if (!(node_info.in_block(position)))
      node_iter.next_skipBranch();
    else
      {
      if ((!(node_iter.children_present())) 
	  || (node_iter.get_nodeID().get_level() >= max_level))
	found_it = true;
      else
	node_iter.next();
      }
    }

  if (!found_it)
    throw Except("Block not found", __FILE__, __LINE__);

  const unsigned int block_index = node_iter.get_data().get_block_index();

  tree_iter = node_ptrs[block_index];
  }


const Tree::NodeData & MeshInfo::get_node_data(const std::vector<double> 
					       & position) const
  {
  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  // Identify the tree for this position

  const unsigned int num_trees = tree_block.get_num_cells();

  bool found_tree = false;

  unsigned int tree_index = 0;

  while ((!found_tree) && (tree_index < num_trees))
    {
    if (tree_block[tree_index].get_data().get_block_info().in_block(position))
      found_tree = true;
    else
      tree_index++;
    }

  if (!found_tree)
    throw Except("Tree not found", __FILE__, __LINE__);

  Tree::ConstTreeNodeIter node_iter = tree_block[tree_index].get_node_iter();

  bool found_it = false;

  while ((!found_it) && (node_iter.valid_node()))
    {
    const Block::BlockInfo & node_info = node_iter.get_data().get_block_info();

    if (!(node_info.in_block(position)))
      node_iter.next_skipBranch();
    else
      {
      if (node_iter.children_present())
        node_iter.next();
      else
        found_it = true;
      }
    }

  if (!found_it)
    throw Except("Block not found", __FILE__, __LINE__);

  return node_iter.get_data();
  }


const Tree::NodeData & MeshInfo::get_node_data(const std::vector<double> 
                                                 & position,
					       const unsigned int max_level) 
  const
  {
  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  // Identify the tree for this position

  const unsigned int num_trees = tree_block.get_num_cells();

  bool found_tree = false;

  unsigned int tree_index = 0;

  while ((!found_tree) && (tree_index < num_trees))
    {
    if (tree_block[tree_index].get_data().get_block_info().in_block(position))
      found_tree = true;
    else
      tree_index++;
    }

  if (!found_tree)
    throw Except("Tree not found", __FILE__, __LINE__);

  Tree::ConstTreeNodeIter node_iter = tree_block[tree_index].get_node_iter();

  bool found_it = false;

  while ((!found_it) && (node_iter.valid_node()))
    {
    const Block::BlockInfo & node_info = node_iter.get_data().get_block_info();

    if (!(node_info.in_block(position)))
      node_iter.next_skipBranch();
    else
      {
      if ((!(node_iter.children_present())) 
	  || (node_iter.get_nodeID().get_level() >= max_level))
	found_it = true;
      else
	node_iter.next();
      }
    }

  if (!found_it)
    throw Except("Block not found", __FILE__, __LINE__);

  return node_iter.get_data();
  }


void
MeshInfo::get_cell_index (
    const std::vector<double> & position,
	unsigned int & block_index,
    unsigned int & cell_index) const
{
    if (!valid_mesh)
        throw Except("No valid mesh defined", __FILE__, __LINE__);

    const Tree::NodeData & node_data = get_node_data(position);

    block_index = node_data.get_block_index();

    const Block::BlockInfo & block_info = node_data.get_block_info();

    cell_index = block_info.get_nearest_cell(position);
}

void
MeshInfo::get_cell_index (
    const std::vector<double> & position,
	unsigned int & block_index,
    unsigned int & cell_index) const
{
    if (!valid_mesh)
        throw Except("No valid mesh defined", __FILE__, __LINE__);

    const Tree::NodeData & node_data = get_node_data(position);

    block_index = node_data.get_block_index();

    const Block::BlockInfo & block_info = node_data.get_block_info();

    cell_index = block_info.get_nearest_cell(position);
}

void MeshInfo::get_cell_index(const std::vector<double> & position,
			      const unsigned int max_level,
			      unsigned int & block_index,
			      unsigned int & cell_index) const
  {
  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  const Tree::NodeData & node_data = get_node_data(position, max_level);

  block_index = node_data.get_block_index();

  const Block::BlockInfo & block_info = node_data.get_block_info();

  cell_index = block_info.get_nearest_cell(position);
  }


void MeshInfo::get_sim_params(const SimInfo & sim_info)
  {
  // Get information from the simulation parameters

  valid_mesh = sim_info.mesh_present();

  dims = sim_info.get_mesh_dims();

  num_blocks = sim_info.get_num_blocks();

  if (num_blocks < 1)
    throw Except("No blocks present in file", __FILE__, __LINE__);

  sim_time = sim_info.get_sim_time();

  block_dims = sim_info.get_block_dims();
  num_block_cells = sim_info.get_cells_per_block();

  base_block_dims = sim_info.get_base_block_dims();

  volume_minbounds = sim_info.get_volume_minbounds();
  volume_maxbounds = sim_info.get_volume_maxbounds();

  volume_center = sim_info.get_volume_center();
  volume_width = sim_info.get_volume_width();

  geom_type = sim_info.get_mesh_geometry();
  }


void MeshInfo::construct_tree(const hid_t file_id, const SimInfo & sim_info)
  {
  // WARNING: Call this function only after calling read_sim_params!

  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  // Get information about the mesh topology

  const FlashVersion flash_version = sim_info.get_flash_code_version();

  const TreeInfo tree_info(file_id, flash_version, dims);

  if (tree_info.get_num_blocks() != num_blocks)
    throw Except("Block count disagreement with SimInfo data", __FILE__, 
		 __LINE__);

  block_refine_level = tree_info.get_block_refine_level();
  block_process_id = tree_info.get_process_id();

  parent_indexes = tree_info.get_block_parent_index();

  const std::vector< std::vector<double> > & block_centers 
    = tree_info.get_block_center();

  const std::vector< std::vector<unsigned int> > & child_indexes
    = tree_info.get_block_child_indexes();

  Tree::build_tree_block(dims, geom_type, volume_minbounds, volume_maxbounds,
			 base_block_dims, block_dims, parent_indexes,
			 block_centers, child_indexes, tree_block);

  neighbor_indexes_lo = tree_info.get_block_neighbor_indexes_lo();
  neighbor_indexes_hi = tree_info.get_block_neighbor_indexes_hi();

  neighbor_types_lo = tree_info.get_block_neighbor_types_lo();
  neighbor_types_hi = tree_info.get_block_neighbor_types_hi();

  boundary_types_lo = tree_info.get_block_boundary_types_lo();

  boundary_types_hi = tree_info.get_block_boundary_types_hi();

  set_node_data_pointers();

  const std::vector<unsigned int> & block_type_ids 
    = tree_info.get_block_type();

  Block::convert_block_type(block_type_ids, block_type);

  // Count the number of leaves and set min/max leaf refine levels

  num_leaf_blocks = 0;

  min_leaf_refine_level = 0;
  max_leaf_refine_level = 0;

  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    if (block_type[block_index] == Block::Leaf)
      {
      // Test the refinement before incrementing num_leaf_blocks

      const unsigned int refine_level = block_refine_level[block_index];

      if (num_leaf_blocks < 1)
	{
	min_leaf_refine_level = refine_level;
	max_leaf_refine_level = refine_level;
	}
      else
	{
	if (refine_level < min_leaf_refine_level)
	  min_leaf_refine_level = refine_level;

	if (refine_level > max_leaf_refine_level)
	  max_leaf_refine_level = refine_level;
	}

      // Count this block

      num_leaf_blocks++;
      }
  }


void MeshInfo::set_node_data_pointers()
  {
  //  Call build_tree_block before calling this function

  if (!valid_mesh)
    throw Except("No valid mesh defined", __FILE__, __LINE__);

  // Initialize pointers to null

  node_ptrs.resize(num_blocks);

  for (unsigned int i = 0; i < num_blocks; i++)
    node_ptrs[i].reset();

  // Traverse the tree

  Tree::ConstTreeIter tree_iter(tree_block);

  while (tree_iter.valid_node())
    {
    const Tree::NodeData & node_data = tree_iter.get_node_data();

    if (!node_data.initialized())
      throw Except("Node not initialized", __FILE__, __LINE__);

    const unsigned int block_index = node_data.get_block_index();

    if (node_ptrs[block_index].valid_node())
      throw Except("Block pointer already assigned", __FILE__, __LINE__);

    node_ptrs[block_index] = tree_iter;

    tree_iter.next();
    }

  if (check_node_data_pointers() < 0)
    throw Except("Unable to assign all node data pointers", __FILE__, 
		 __LINE__);
  }


int MeshInfo::check_node_data_pointers() const
  {
  int ret_val = 0;  // OK

  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    if (!(node_ptrs[block_index].valid_node()))
      ret_val--;  // Decrement below zero to accumulate number missing

  return ret_val;
  }


}  // End namespace File
}  // End namespace QuickFlash
