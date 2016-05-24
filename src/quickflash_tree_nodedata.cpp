// C++ program file quickflash_tree_nodedata.cpp

/*
  By Nathan C. Hearn
     August 10, 2007

  Block tree node data class.
*/


#include "quickflash_tree_nodedata.hpp"
#include <vector>
#include "quickflash_geometry.hpp"


namespace QuickFlash
{
namespace Tree
{

// Class NodeData

NodeData::NodeData() : 
  initialized_data(false), blockindex(0), blockinfo()
  { reset(); }


NodeData::NodeData(const unsigned int block_index, 
		   const std::vector<double> & block_mincoords,
		   const std::vector<double> & block_maxcoords,
		   const std::vector<unsigned int> & block_dims,
		   const Geometry::GeometryType geom_type) :
  initialized_data(false), blockindex(0), blockinfo()
  {
  reset(block_index, block_mincoords, block_maxcoords, block_dims,
	geom_type);
  }


NodeData::NodeData(const NodeData & source) :
  initialized_data(false), blockindex(0), blockinfo()
  { reset(source); }


void NodeData::reset()
  {
  blockindex = 0;

  blockinfo.reset();

  initialized_data = false;
  }


void NodeData::reset(const unsigned int block_index, 
		     const std::vector<double> & block_mincoords,
		     const std::vector<double> & block_maxcoords,
		     const std::vector<unsigned int> & block_dims, 
		     const Geometry::GeometryType geom_type)
  {
  blockindex = block_index;

  blockinfo.reset(block_mincoords, block_maxcoords, block_dims, geom_type);

  initialized_data = true;
  }


void NodeData::reset(const NodeData & source)
  {
  if (&source != this)
    {
    initialized_data = source.initialized_data;

    blockindex = source.blockindex;
    
    blockinfo.reset(source.blockinfo);
    }
  }


}  // End namespace Tree
}  // End namespace QuickFlash
