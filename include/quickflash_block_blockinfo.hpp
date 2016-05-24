// C++ header file quickflash_block_blockinfo.hpp

/*
  By Nathan C. Hearn
     September 2, 2006

  Meta-info for block.
*/


#ifndef QUICKFLASH_BLOCK_BLOCKINFO_HPP
#define QUICKFLASH_BLOCK_BLOCKINFO_HPP


#include <vector>
#include "quickflash_block_utils.hpp"
#include "quickflash_geometry.hpp"
#include "quickflash_block_cellvolume.hpp"


namespace QuickFlash
{
namespace Block
{

//! Block meta-information (spatial coordinates, widths, etc.) class
class BlockInfo
  {
  public:

    //! Default constructor
    BlockInfo() :
      dims(0), numcells(0), blockdims(), mincoords(), maxcoords(), 
      cellwidth(), oneOver_cellwidth(),
      cellvolume(0.0), geom_type(Geometry::Cartesian), phys_cell_vol(0)
      { reset(); }

    //! Initialization with C-style arrays
    BlockInfo(const unsigned int space_dims, const double block_mincoords[],
	      const double block_maxcoords[], 
	      const unsigned int block_dims[],
	      const Geometry::GeometryType geometry_type) :
      dims(0), numcells(0), blockdims(), mincoords(), maxcoords(), 
      cellwidth(), oneOver_cellwidth(),
      cellvolume(0.0), geom_type(geometry_type), phys_cell_vol(0)
      { 
      reset(space_dims, block_mincoords, block_maxcoords, block_dims,
	    geometry_type); 
      }

    //! Initialization with STL vectors
    BlockInfo(const std::vector<double> & block_mincoords,
	      const std::vector<double> & block_maxcoords,
	      const std::vector<unsigned int> & block_dims,
	      const Geometry::GeometryType geometry_type) :
      dims(0), numcells(0), blockdims(), mincoords(), maxcoords(), 
      cellwidth(), oneOver_cellwidth(),
      cellvolume(0.0), geom_type(geometry_type), phys_cell_vol(0)
      { reset(block_mincoords, block_maxcoords, block_dims, geometry_type); }

    //! Copy constructor
    BlockInfo(const BlockInfo & source) :
      dims(0), numcells(0), blockdims(), mincoords(), maxcoords(), 
      cellwidth(), oneOver_cellwidth(),
      cellvolume(0.0), geom_type(Geometry::Cartesian), phys_cell_vol(0)
      { reset(source); }

    //! Destructor
    ~BlockInfo() 
      { 
      reset(); // Deletes phys_cell_vol 
      }

    //! Assignment operator
    BlockInfo & operator=(const BlockInfo & source)
      {
      reset(source);
      return *this;
      }

    //! Reset block to default (uninitialized) values
    void reset();

    //! Reset block information with C-style arrays
    void reset(const unsigned int space_dims, const double block_mincoords[],
	       const double block_maxcoords[], 
	       const unsigned int block_dims[],
	       const Geometry::GeometryType geometry_type);

    //! Reset block information with STL vectors
    void reset(const std::vector<double> & block_mincoords,
	       const std::vector<double> & block_maxcoords,
	       const std::vector<unsigned int> & block_dims,
	       const Geometry::GeometryType geometry_type);

    //! Copy block information from existing BlockInfo object
    void reset(const BlockInfo & source);

    //! Returns the number of spatial dimensions
    unsigned int get_space_dims() const { return dims; }

    //! Returns the number of cells along each axis
    const std::vector<unsigned int> & get_block_dims() const 
      { return blockdims; }

    //! Returns the minimum coordinates of the block
    const std::vector<double> & get_min_coords() const { return mincoords; }

    //! Returns the maximum coordinates of the block
    const std::vector<double> & get_max_coords() const { return maxcoords; }

    //! Computes the center coordinates of the block
    void get_block_center(std::vector<double> & block_center) const 
      {
      block_center.resize(dims);

      for (unsigned int i = 0; i < dims; i++)
	block_center[i] = 0.5 * (mincoords[i] + maxcoords[i]);
      }

    //! Computes the width of the block
    void get_block_width(std::vector<double> & block_width) const 
      {
      block_width.resize(dims);

      for (unsigned int i = 0; i < dims; i++)
	block_width[i] = maxcoords[i] - mincoords[i];
      }

    //! Returns the width of a cell along each axis
    const std::vector<double> & get_cell_width() const { return cellwidth; }

    //! Returns the reciprical of the width of a cell along each axis
    const std::vector<double> & get_oneOver_cell_width() const 
      { return oneOver_cellwidth; }

    //! Returns the volume of a cell (no geometry information
    double get_cell_volume() const { return cellvolume; }

    //! Returns the true (geometric) volume of a cell
    double get_cell_volume(const unsigned int cell_index) const
      { return phys_cell_vol->cell_volume(cell_index); }

    //! Returns the true (geometric) volume of a cell
    double get_cell_volume(const std::vector<unsigned int> & cell_index) const
      { return phys_cell_vol->cell_volume(cell_index); }

    //! Returns the number of cells in the block
    unsigned int get_num_cells() const { return numcells; }

    //! Returns the mesh geometry for the block
    Geometry::GeometryType get_geometry_type() const { return geom_type; }

    //! Returns true if coorindates are less than maxcoords but not mincoords
    bool in_block(const std::vector<double> & coords) const
      {
      bool inside = true;

      for (unsigned int i = 0; i < dims; i++)
	{
	const double axis_coord = coords[i];

	if ((axis_coord < mincoords[i]) || (!(axis_coord < maxcoords[i])))
	  inside = false;
	}

      return inside;
      }

    //! Returns true if cell_index is a valid coordinate in this lbock
    bool in_block(const std::vector<unsigned int> & cell_index) const
      {
      bool inside = true;

      for (unsigned int i = 0; i < dims; i++)
	{
	const double axis_index = cell_index[i];

	if ((axis_index < 0) || (!(axis_index < blockdims[i])))
	  inside = false;
	}

      return inside;
      }


    //! Returns true if coorindates are within the cell (includes min boundary)
    bool in_cell(const unsigned int cell_index,
		 const std::vector<double> & coords) const
      {
      std::vector<double> cell_min(dims);
      std::vector<double> cell_max(dims);

      get_cell_bounds(cell_index, cell_min, cell_max);

      bool inside = true;

      for (unsigned int i = 0; i < dims; i++)
	{
	const double axis_coord = coords[i];

	if ((axis_coord < cell_min[i]) || (!(axis_coord < cell_max[i])))
	  inside = false;
	}

      return inside;
      }


    //! Returns true if coorindates are within the cell (includes min boundary)
    bool in_cell(const std::vector<unsigned int> & cell_index,
		 const std::vector<double> & coords) const
      {
      std::vector<double> cell_min(dims);
      std::vector<double> cell_max(dims);

      get_cell_bounds(cell_index, cell_min, cell_max);

      bool inside = true;

      for (unsigned int i = 0; i < dims; i++)
	{
	const double axis_coord = coords[i];

	if ((axis_coord < cell_min[i]) || (!(axis_coord < cell_max[i])))
	  inside = false;
	}

      return inside;
      }


    //! Returns the coordinates of the center of the specified cell (scalar)
    void get_cell_center(const unsigned int cell_index,
			std::vector<double> & center_coords) const
      {
      std::vector<unsigned int> vector_index(dims);
      get_vector_index(cell_index, vector_index);

      get_cell_center(vector_index, center_coords);
      }

    //! Returns the coordinates of the center of the specified cell (vector)
    void get_cell_center(const std::vector<unsigned int> & cell_index,
			 std::vector<double> & center_coords) const
      {
      center_coords = mincoords;

      for (unsigned int i = 0; i < dims; i++)
	center_coords[i] += cellwidth[i] * 
	  (static_cast<double>(cell_index[i]) + 0.5);
      }

    //! Returns the minimum and maximum bounds of the specified cell (scalar)
    void get_cell_bounds(const unsigned int cell_index,
			 std::vector<double> & min_coords,
			 std::vector<double> & max_coords) const
      {
      std::vector<unsigned int> vector_index(dims);
      get_vector_index(cell_index, vector_index);

      get_cell_bounds(vector_index, min_coords, max_coords);
      }

    //! Returns the minimum and maximum bounds of the specified cell (vector)
    void get_cell_bounds(const std::vector<unsigned int> & cell_index,
			 std::vector<double> & min_coords,
			 std::vector<double> & max_coords) const;

    //! Returns the index of the cell with a center closest to coords
    unsigned int get_nearest_cell(const std::vector<double> & coords) const;

    //! Returns the index of the cell with a center closest to coords
    void get_nearest_cell(const std::vector<double> & coords,
			  std::vector<unsigned int> & cell_index) const;

    //! Converts from a scalar cell index to a vector cell index
    void get_vector_index(const unsigned int scalar_index,
			  std::vector<unsigned int> & vector_index) const;

    //! Converts from a vector cell index to a scalar cell index
    unsigned int get_scalar_index(const std::vector<unsigned int> 
				  & vector_index) const;

    //! Inverts the cell index (reflects across the midpoint) along an axis
    void invert_cell_index(const unsigned int axis, 
			   const std::vector<unsigned int> & cell_index, 
			   std::vector<unsigned int> & inverted_cell_index)
      const
      { 
      Utils::invert_cell_index(blockdims, axis, cell_index, 
			       inverted_cell_index);
      }

    //! Inverts the cell index (reflects across the midpoint) along an axis,
    //! over-writing the index vector
    void invert_cell_index_in_place(const unsigned int axis, 
				    std::vector<unsigned int> & cell_index)
      const
      { Utils::invert_cell_index_in_place(blockdims, axis, cell_index); }

    //! Inverts the cell index (reflects across the midpoint) along an axis
    unsigned int invert_cell_index(const unsigned int axis, 
				   const unsigned int cell_index) const
      { return Utils::invert_cell_index(blockdims, axis, cell_index); }

    //! Determine the boundaries of a specific child oct-tree node
    void get_child_bounds(const unsigned int local_child_index,
			  std::vector<double> & child_minbounds,
			  std::vector<double> & child_maxbounds) const
      {
      Utils::get_subblock_bounds(local_child_index, mincoords, maxcoords, 
				 child_minbounds, child_maxbounds);
      }

    unsigned int get_nearest_cell_axis_index(const unsigned int axis, 
					     const double axis_pos) const;
		     
  private:

    void compute_cell_volume();

    void get_axis_cell_bounds(const unsigned int axis,
			      const unsigned int axis_cell_index,
			      double & min_bound, double & max_bound) const;

  private:

    unsigned int dims;
    unsigned int numcells;
    std::vector<unsigned int> blockdims;

    std::vector<double> mincoords;
    std::vector<double> maxcoords;

    std::vector<double> cellwidth;
    std::vector<double> oneOver_cellwidth;
    double cellvolume;

    Geometry::GeometryType geom_type;

    const CellVolume::VolumeBase * phys_cell_vol;
  };


}  // End namespace Block
}  // End namespace QuickFlash


#endif  // QUICKFLASH_BLOCK_BLOCKINFO_HPP
