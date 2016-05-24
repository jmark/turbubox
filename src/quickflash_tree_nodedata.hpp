// C++ header file quickflash_tree_nodedata.hpp

/*
  By Nathan C. Hearn
     October 27, 2006

  Block tree node data class.
*/


#ifndef QUICKFLASH_TREE_NODEDATA_HPP
#define QUICKFLASH_TREE_NODEDATA_HPP


#include <vector>
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_geometry.hpp"


namespace QuickFlash
{
namespace Tree
{

// Class NodeData

class NodeData
  {
  public:
    NodeData();

    NodeData(const unsigned int block_index, 
	     const std::vector<double> & block_mincoords,
	     const std::vector<double> & block_maxcoords,
	     const std::vector<unsigned int> & block_dims,
	     const Geometry::GeometryType geom_type);

    NodeData(const NodeData & source);

    ~NodeData() { }

    NodeData & operator=(const NodeData & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const unsigned int block_index, 
	       const std::vector<double> & block_mincoords,
	       const std::vector<double> & block_maxcoords,
	       const std::vector<unsigned int> & block_dims, 
	       const Geometry::GeometryType geom_type);

    void reset(const NodeData & source);

    bool initialized() const { return initialized_data; }

    unsigned int get_block_index() const { return blockindex; }

    const Block::BlockInfo & get_block_info() const { return blockinfo; }

    void get_child_bounds(const unsigned int local_child_index,
			  std::vector<double> & child_minbounds,
			  std::vector<double> & child_maxbounds) const
    { 
    blockinfo.get_child_bounds(local_child_index, child_minbounds, 
			       child_maxbounds);
    }

  private:
    bool initialized_data;

    unsigned int blockindex;

    Block::BlockInfo blockinfo;
  };


}  // End namespace Tree
}  // End namespace QuickFlash


#endif  // QUICKFLASH_TREE_NODEDATA_HPP
