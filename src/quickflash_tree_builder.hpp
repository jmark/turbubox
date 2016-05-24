// C++ header file quickflash_tree_builder.hpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Tree construction utilities.
*/


#ifndef QUICKFLASH_TREE_BUILDER
#define QUICKFLASH_TREE_BUILDER


#include <vector>
#include "quickflash_tree_block.hpp"
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
		TreeNode & root_node);

void reset_tree(const TreeNode & source_root, TreeNode & receiver_root);

void build_tree_block(const unsigned int dims, 
		      const Geometry::GeometryType geom_type,
		      const std::vector<double> & volume_minbounds, 
		      const std::vector<double> & volume_maxbounds,
		      const std::vector<unsigned int> & base_block_dims,
		      const std::vector<unsigned int> & block_dims,
		      const std::vector<unsigned int> & parent_indexes,
		      const std::vector< std::vector<double> > & block_centers,
	      const std::vector< std::vector<unsigned int> > & child_indexes,
		      TreeBlock & tree_block);

void reset_tree_block(const TreeBlock & source_block, 
		      TreeBlock & receiver_block);


}  // End namespace Tree
}  // End namespace QuickFlash


#endif  // QUICKFLASH_TREE_BUILDER
