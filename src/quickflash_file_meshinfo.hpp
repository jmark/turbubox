// C++ header file quickflash_file_meshinfo.hpp

/*
  By Nathan C. Hearn
     October 4, 2008
*/


#ifndef QUICKFLASH_FILE_MESHINFO_HPP
#define QUICKFLASH_FILE_MESHINFO_HPP

#include <hdf5.h>
#include <vector>
#include <string>
#include "quickflash_file_siminfo.hpp"
#include "quickflash_block_utils.hpp"
#include "quickflash_block_blocktype.hpp"
#include "quickflash_tree_block.hpp"
#include "quickflash_tree_nodedata.hpp"
#include "quickflash_geometry.hpp"
#include "quickflash_mesh_boundarydefs.hpp"
#include "quickflash_file_flashdefs.hpp"


namespace QuickFlash
{
namespace File
{

// Class MeshInfo

class MeshInfo
  {
  public :

    MeshInfo();

    MeshInfo(const hid_t file_id, const SimInfo & sim_info);

    MeshInfo(const MeshInfo & source);

    ~MeshInfo() { }

    MeshInfo & operator=(const MeshInfo & source)
      {
      reset(source);
      return *this;
      }

    void reset();
    void reset(const hid_t file_id, const SimInfo & sim_info);
    void reset(const MeshInfo & source);

    bool mesh_info_read() const { return mesh_info_set_up; }

    bool mesh_valid() const { return valid_mesh; }

    double get_sim_time() const { return sim_time; }

    unsigned int get_dims() const { return dims; }

    unsigned int get_num_blocks() const { return num_blocks; }

    unsigned int get_num_leaf_blocks() const { return num_leaf_blocks; }

    const std::vector<unsigned int> & get_base_block_dims() const
      { return base_block_dims; }

    const std::vector<double> & get_volume_minbounds() const
      { return volume_minbounds; }

    const std::vector<double> & get_volume_maxbounds() const
      { return volume_maxbounds; }

    const std::vector<double> & get_volume_center() const
      { return volume_center; }

    const std::vector<double> & get_volume_width() const
      { return volume_width; }

    bool in_bounds(const std::vector<double> & position) const;

    bool is_leaf(const unsigned int block_index) const
      { return (block_type[block_index] == Block::Leaf); }

    bool is_root(const unsigned int block_index) const
      { return (get_parent_index(block_index) == block_index); }

    Tree::ConstTreeIter get_tree_iter() const 
      { return Tree::ConstTreeIter(tree_block); }

    void get_tree_iter(Tree::ConstTreeIter & tree_iter) const
      { tree_iter.reset(tree_block); }

    Tree::ConstTreeIter get_tree_iter(const unsigned int block_index) const
      { return Tree::ConstTreeIter(node_ptrs[block_index]); }

    void get_tree_iter(const unsigned int block_index, 
		       Tree::ConstTreeIter & tree_iter) const
      { tree_iter.reset(node_ptrs[block_index]); }

    Tree::ConstTreeIter get_tree_iter(const std::vector<double> & position) 
      const
      {
      Tree::ConstTreeIter tree_iter;

      get_tree_iter(position, tree_iter);

      return Tree::ConstTreeIter(tree_iter);
      }

    Tree::ConstTreeIter get_tree_iter(const std::vector<double> & position,
				      const unsigned int max_level) 
      const
      {
      Tree::ConstTreeIter tree_iter;

      get_tree_iter(position, max_level, tree_iter);

      return Tree::ConstTreeIter(tree_iter);
      }

    void get_tree_iter(const std::vector<double> & position,
		       Tree::ConstTreeIter & tree_iter) const;

    void get_tree_iter(const std::vector<double> & position,
		       const unsigned int max_level,
		       Tree::ConstTreeIter & tree_iter) const;

    const Tree::NodeData & get_node_data(const std::vector<double>
					 & position) const;

    const Tree::NodeData & get_node_data(const std::vector<double>
					 & position,
					 const unsigned int max_level) const;

    unsigned int get_block_index(const std::vector<double> & position) const
      { return get_node_data(position).get_block_index(); }

    unsigned int get_block_index(const std::vector<double> & position,
				 const unsigned int max_level) const
      { return get_node_data(position, max_level).get_block_index(); }

    const Block::BlockInfo & get_block_info(const std::vector<double>
					    & position) const
      { return get_node_data(position).get_block_info(); }

    const Block::BlockInfo & get_block_info(const std::vector<double>
					    & position,
					    const unsigned int max_level) const
      { return get_node_data(position, max_level).get_block_info(); }

    void get_cell_index(const std::vector<double> & position,
			unsigned int & block_index,
			unsigned int & cell_index) const;

    void get_cell_index(const std::vector<double> & position,
			const unsigned int max_level,
			unsigned int & block_index,
			unsigned int & cell_index) const;

    Block::BlockType get_block_type(const unsigned int block_index) const
      { return block_type[block_index]; }

    unsigned int get_refine_level(const unsigned int block_index) const
      { return block_refine_level[block_index]; }

    unsigned int get_min_leaf_refine_level() const 
      { return min_leaf_refine_level; }

    unsigned int get_max_leaf_refine_level() const 
      { return max_leaf_refine_level; }

    unsigned int get_process_id(const unsigned int block_index) const
      { return block_process_id[block_index]; }

    const Tree::NodeData & get_node_data(const unsigned int block_index) const
      { return node_ptrs[block_index].get_node_data(); }

    const Block::BlockInfo & get_block_info(const unsigned int block_index)
      const
      { return get_node_data(block_index).get_block_info(); }

    void get_block_center(const unsigned int block_index, 
			  std::vector<double> & block_center) const
      { get_block_info(block_index).get_block_center(block_center); }

    void get_block_width(const unsigned int block_index, 
			 std::vector<double> & block_width) const
      { get_block_info(block_index).get_block_width(block_width); }

    const std::vector<double> & get_block_mincoords(const unsigned int 
						    block_index) const
      { return get_block_info(block_index).get_min_coords(); }

    const std::vector<double> & get_block_maxcoords(const unsigned int 
						    block_index) const
      { return get_block_info(block_index).get_max_coords(); }

    const std::vector<unsigned int> & get_block_dims() const 
      { return block_dims; }

    unsigned int get_num_block_cells() const { return num_block_cells; }

    unsigned int get_parent_index(const unsigned int block_index) const
      { return parent_indexes[block_index]; }

    const std::vector<unsigned int> & 
      get_neighbor_indexes_lo(const unsigned int block_index) const
      { return neighbor_indexes_lo[block_index]; }

    const std::vector<unsigned int> & 
      get_neighbor_indexes_hi(const unsigned int block_index) const
      { return neighbor_indexes_hi[block_index]; }

    const std::vector<Mesh::NeighborType> & 
      get_neighbor_types_lo(const unsigned int block_index) const
      { return neighbor_types_lo[block_index]; }

    const std::vector<Mesh::NeighborType> & 
      get_neighbor_types_hi(const unsigned int block_index) const
      { return neighbor_types_hi[block_index]; }

    const std::vector<Mesh::BoundaryType> & 
      get_boundary_types_lo(const unsigned int block_index) const
      { return boundary_types_lo[block_index]; }

    const std::vector<Mesh::BoundaryType> & 
      get_boundary_types_hi(const unsigned int block_index) const
      { return boundary_types_hi[block_index]; }

    void get_sub_block_bounds(const unsigned int block_index,
			      const unsigned int child_index,
			      std::vector<double> & min_bounds,
			      std::vector<double> & max_bounds)
      { 
      Block::Utils::get_subblock_bounds(child_index, 
					get_block_mincoords(block_index),
					get_block_maxcoords(block_index), 
					min_bounds, max_bounds);
      }

  private :

    void get_sim_params(const SimInfo & sim_info);

    void construct_tree(const hid_t file_id, const SimInfo & sim_info);

    void set_node_data_pointers();
    int check_node_data_pointers() const;

  private :

    bool mesh_info_set_up;

    bool valid_mesh;

    double sim_time;

    unsigned int dims;

    unsigned int num_blocks;
    unsigned int num_leaf_blocks;

    unsigned int min_leaf_refine_level;
    unsigned int max_leaf_refine_level;

    std::vector<unsigned int> block_refine_level;
    std::vector<unsigned int> block_process_id;

    std::vector<unsigned int> base_block_dims;

    std::vector<double> volume_minbounds;
    std::vector<double> volume_maxbounds;

    std::vector<double> volume_center;
    std::vector<double> volume_width;

    std::vector<unsigned int> block_dims;
    unsigned int num_block_cells;

    Tree::TreeBlock tree_block;

    std::vector<Tree::ConstTreeIter> node_ptrs;

    std::vector<Block::BlockType> block_type;

    Geometry::GeometryType geom_type;

    std::vector<unsigned int> parent_indexes;

    std::vector< std::vector<unsigned int> > neighbor_indexes_lo;
    std::vector< std::vector<unsigned int> > neighbor_indexes_hi;

    std::vector< std::vector<Mesh::NeighborType> > neighbor_types_lo;
    std::vector< std::vector<Mesh::NeighborType> > neighbor_types_hi;

    std::vector< std::vector<Mesh::BoundaryType> > boundary_types_lo;
    std::vector< std::vector<Mesh::BoundaryType> > boundary_types_hi;
  };


}  // End namespace File
}  // End namespace QuickFlash


#endif  // QUICKFLASH_FILE_MESHINFO_HPP
