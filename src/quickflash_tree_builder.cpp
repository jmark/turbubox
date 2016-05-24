// C++ program file quickflash_tree_builder.cpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Tree construction utilities.
*/


#include "quickflash_tree_builder.hpp"
#include <vector>
#include "quickflash_tree_nodedata.hpp"
#include "quickflash_except.hpp"
#include "quickflash_tree_block.hpp"
#include "quickflash_tree_nodeid.hpp"
#include "quickflash_geometry.hpp"


namespace QuickFlash
{
namespace Tree
{

void build_tree(const unsigned int dims, 
		const Geometry::GeometryType geom_type,
		const std::vector<unsigned int> & block_dims, 
		const std::vector< std::vector<double> > & block_centers,
		const std::vector< std::vector<unsigned int> > & child_indexes,
		TreeNode & root_node)
  {
  /*
    Assume two things:

      1. The root node for this tree has already been initialized
      2. block_info has the correct length (num_blocks)
  */

  TreeNode * node_ptr = root_node.get_nodePtr();

  std::vector<double> temp_minbounds(dims);
  std::vector<double> temp_maxbounds(dims);

  while (node_ptr != 0)
    {
    NodeData & node_data = node_ptr->get_data();

    if (!(node_data.initialized()))
      throw Except("Tree node not initialized", __FILE__, __LINE__);

    const unsigned int current_block_index = node_data.get_block_index();

    const std::vector<double> & current_block_center
      = block_centers[current_block_index];

    // Store the node's child block indexes

    const std::vector<unsigned int> & child_index_list
      = child_indexes[current_block_index];

    // Create the children

    const unsigned int num_children = child_index_list.size();

    for (unsigned int child = 0; child < num_children; child++)
      {
      const unsigned int child_index = child_index_list[child];

      // Get the center coordinates for the child
  
      const std::vector<double> & child_center = block_centers[child_index];

      // Determine the scalar child-block index for this node

      unsigned int scalar_index = 0;

      for (unsigned int i = 0; i < dims; i++)
	if (!(child_center[i] < current_block_center[i]))
	  scalar_index += 1 << i;

      // Create this child

      if (node_ptr->child_present(scalar_index))
	throw Except("Coincident child nodes", __FILE__, __LINE__);

      TreeNode * child_ptr = node_ptr->create_child(scalar_index);

      NodeData & child_data = child_ptr->get_data();

      node_data.get_child_bounds(scalar_index, temp_minbounds, 
				 temp_maxbounds);

      child_data.reset(child_index, temp_minbounds, temp_maxbounds,
		       block_dims, geom_type);
      }

    node_ptr = node_ptr->get_nextNode();
    }  
  }


void reset_tree(const TreeNode & source_root, TreeNode & receiver_root)
  {
  receiver_root.reset();

  const TreeNode * source_node = source_root.get_nodePtr();

  TreeNode * recv_node = receiver_root.get_nodePtr();

  while ((source_node != 0) && (recv_node != 0))
    {
    // ANY ERROR CHECKING NEEDED???

    // Copy the data

    const NodeData & source_data = source_node->get_data();

    recv_node->get_data().reset(source_data);

    // Create the children

    const TreeNode * source_child = source_node->get_firstChild();

    while (source_child != 0)
      {
      const NodeID & child_id = source_child->get_nodeID();

      const LevelID_Type child_level_id = child_id.get_level_id();

      recv_node->create_child(child_level_id);

      source_child = source_child->get_nextSibling();
      }

    source_node = source_node->get_nextNode();
    recv_node = recv_node->get_nextNode();
    }
  }


void build_tree_block(const unsigned int dims, 
		      const Geometry::GeometryType geom_type,
		      const std::vector<double> & volume_minbounds, 
		      const std::vector<double> & volume_maxbounds,
		      const std::vector<unsigned int> & base_block_dims,
		      const std::vector<unsigned int> & block_dims,
		      const std::vector<unsigned int> & parent_indexes,
		      const std::vector< std::vector<double> > & block_centers,
	      const std::vector< std::vector<unsigned int> > & child_indexes,
		      TreeBlock & tree_block)
  {
  // Assume that block_info already has the correct length (num_blocks)

  // Set up the root nodes

  tree_block.reset(volume_minbounds, volume_maxbounds, base_block_dims,
		   geom_type);
  
  const unsigned int num_blocks = block_centers.size();

  for (unsigned int index = 0; index < num_blocks; index++)
    {
    // A root node is one whose parent is itself

    if (parent_indexes[index] == index)
      {
      const std::vector<double> & cell_center = block_centers[index];

      const unsigned int cell_index = tree_block.get_nearest_cell(cell_center);

      tree_block[cell_index].reset();

      std::vector<double> cell_mincoords;
      std::vector<double> cell_maxcoords;

      tree_block.get_cell_bounds(cell_index, cell_mincoords, cell_maxcoords);

      NodeData & node_data = tree_block[cell_index].get_data();

      if (node_data.initialized())
	throw Except("Overlapping root nodes", __FILE__, __LINE__);

      node_data.reset(index, cell_mincoords, cell_maxcoords, block_dims,
		      geom_type);

      // Build this tree

      build_tree(dims, geom_type, block_dims, block_centers, child_indexes, 
		 tree_block[cell_index]);
      }
    }
  }


void reset_tree_block(const TreeBlock & source_block, 
		      TreeBlock & receiver_block)
  {
  const std::vector<double> & min_coords = source_block.get_min_coords();
  const std::vector<double> & max_coords = source_block.get_max_coords();

  const std::vector<unsigned int> & mesh_dims = source_block.get_block_dims();

  const Geometry::GeometryType geom_type = source_block.get_geometry_type();

  receiver_block.reset(min_coords, max_coords, mesh_dims, geom_type);

  // ANY ERROR CHECKING NEEDED???

  const unsigned int num_cells = source_block.get_num_cells();

  for (unsigned int cell_index = 0; cell_index < num_cells; cell_index++)
    {
    const TreeNode & source_tree = source_block[cell_index];
    TreeNode & receiver_tree = receiver_block[cell_index];

    receiver_tree.set_data(source_tree.get_data());

    reset_tree(source_tree, receiver_tree);
    }
  }


}  // End namespace Tree
}  // End namespace QuickFlash
