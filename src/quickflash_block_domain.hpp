// C++ header file quickflash_block_domain.hpp
/*
  By Nathan C. Hearn
     October 27, 2006

  Template for regularly-spaced volume data storage.

  NOTE: Used in quickflash_tree_block.hpp.
*/


#ifndef QUICKFLASH_BLOCK_DOMAIN_HPP
#define QUICKFLASH_BLOCK_DOMAIN_HPP


#include <vector>
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_geometry.hpp"


namespace QuickFlash
{
namespace Block
{

template<class T>
class Domain : public BlockInfo
  {
  public:
    Domain() : BlockInfo(), mesh() { reset(); }

    Domain(const unsigned int space_dims, const double min_coords[], 
	   const double max_coords[], const unsigned int mesh_dims[],
	   const Geometry::GeometryType geom_type) :
      BlockInfo(), mesh()
      { reset(space_dims, min_coords, max_coords, mesh_dims, geom_type); }

    Domain(const std::vector<double> & min_coords,
	   const std::vector<double> & max_coords,
	   const std::vector<unsigned int> & mesh_dims,
	   const Geometry::GeometryType geom_type) :
      BlockInfo(), mesh()
      { reset(min_coords, max_coords, mesh_dims, geom_type); }

    Domain(const Domain<T> & source) :
      BlockInfo(), mesh()
      { reset(source); }

    ~Domain() { }

    Domain<T> & operator=(const Domain<T> & source)
      {
      reset(source);
      return *this;
      }

    void reset() 
      {
      BlockInfo::reset();
      mesh.clear();
      }

    void reset(const unsigned int space_dims, const double min_coords[], 
	       const double max_coords[], const unsigned int mesh_dims[],
	       const Geometry::GeometryType geom_type)
      {
      BlockInfo::reset(space_dims, min_coords, max_coords, mesh_dims,
		       geom_type);
      mesh.resize(get_num_cells());
      }

    void reset(const std::vector<double> & min_coords,
	       const std::vector<double> & max_coords,
	       const std::vector<unsigned int> & mesh_dims,
	       const Geometry::GeometryType geom_type)
      {
      BlockInfo::reset(min_coords, max_coords, mesh_dims, geom_type);
      mesh.resize(get_num_cells());
      }

    void reset(const Domain<T> & source)
      {
      BlockInfo::reset(source);
      mesh = source.mesh;
      }

    T & operator[](const std::vector<unsigned int> & cell_index)
      { return mesh[get_scalar_index(cell_index)]; }

    const T & operator[](const std::vector<unsigned int> & cell_index) const
      { return mesh[get_scalar_index(cell_index)]; }

    T & operator[](const unsigned int cell_index) 
      { return mesh[cell_index]; }

    const T & operator[](const unsigned int cell_index) const
      { return mesh[cell_index]; }

  protected:
    std::vector<T> mesh;
  };


}  // End namespace Block
}  // End namespace QuickFlash


#endif  // QUICKFLASH_BLOCK_DOMAIN_HPP
