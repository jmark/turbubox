// C++ header file quickflash_block_blockdata.hpp

/*
  By Nathan C. Hearn
     September 2, 2006

  Meta-info for block.
*/


#ifndef QUICKFLASH_BLOCK_BLOCKDATA_HPP
#define QUICKFLASH_BLOCK_BLOCKDATA_HPP


#include <vector>
#include "quickflash_except.hpp"
#include "quickflash_counters.hpp"


namespace QuickFlash
{
namespace Block
{

//! Basic, multi-dimensional block data container template
template<typename T>
class BlockData
  {
  public:
    //! Block data default constructor
    BlockData() : numcells(0), spacedims(0), blockdims(), blockdata() { }

    //! Block constructor for a given spatial dimension and block count
    BlockData(const unsigned int space_dims, const unsigned int block_dims[]) :
      numcells(0), spacedims(0), blockdims(), blockdata()
      { reset(space_dims, block_dims); }

    //! Block constructor fora given spatial dimension and block count (STL)
    BlockData(const std::vector<unsigned int> & block_dims) :
      numcells(0), spacedims(0), blockdims(), blockdata()
      { reset(block_dims); }

    //! Block data copy constructor
    BlockData(const BlockData<T> & source) :
      numcells(0), spacedims(0), 
      blockdims(), blockdata()
      { reset(source); }

    //! Block data destructor
    ~BlockData() { }

    //! Block data assignment operator
    BlockData<T> & operator=(const BlockData<T> & source)
      {
      reset(source);
      return *this;
      }

    //! Returns the number of spatial dimensions
    unsigned int get_space_dims() const { return spacedims; }

    //! Returns the size of the block (in cells) along each axis
    const std::vector<unsigned int> & get_block_dims() const 
      { return blockdims; }

    //! Finds the minimum and maximum values for the block
    void get_minmax_values(T & min_value, T & max_value) const;

    //! Array access operator for cell data (by integer index)
    T & operator[](const unsigned int cell_index) 
      { return blockdata[cell_index]; }

    //! Array access operator for cell data (by integer index) - constant
    const T & operator[](const unsigned int cell_index) const
      { return blockdata[cell_index]; }

    //! Array access operator for cell data (by cell index vector)
    T & operator[](const std::vector<unsigned int> & cell_index)
      { return blockdata[get_scalar_index(cell_index)]; }

    //! Array access operator for cell data (by cell index vector) - constant
    const T & operator[](const std::vector<unsigned int> cell_index) const
      { return blockdata[get_scalar_index(cell_index)]; }

    //! Get internal 1D data array
    const std::vector<T> & get_data() const { return blockdata; }

    //! Replaces current data -- number of cells must match
    void set_data(const std::vector<T> & block_data)
      {
      if (block_data.size() != numcells)
	throw Except("Incompatible data vector", __FILE__, __LINE__);

      blockdata = block_data;
      }

    //! Replaces current data -- number of cells must match
    void set_data(const unsigned int num_elems, const T block_data[])
      {
      if (num_elems != blockdata.size())
	throw Except("Incompatible data vector", __FILE__, __LINE__);

      for (unsigned int i = 0; i < num_elems; i++)
	blockdata[i] = block_data[i];
      }

    //! Replaces current data, converting from Fortran-style rank order
    void set_data_reverse(const unsigned int num_elems, const T block_data[]);

    //! Reset the block dimensions to an un-initialized state
    void reset()
      {
      numcells = 0;
      spacedims = 0;
      blockdims.clear();
      blockdata.clear();
      }

    //! Reset the block dimensions (does not clear data)
    void reset(const unsigned int space_dims, const unsigned int block_dims[]);

    //! Reset the block dimensions (does not clear data)
    void reset(const std::vector<unsigned int> & block_dims);

    //! Copy the dimensions and data from another BlockData object
    void reset(const BlockData<T> & source)
      {
      if (&source != this)
	{
	numcells = source.numcells;
	spacedims = source.spacedims;

	blockdims = source.blockdims;

	blockdata = source.blockdata;
	}
      }

    //! Returns the number of cells in the block
    unsigned int get_num_cells() const { return numcells; }

    //! Convert a vector cell index into a scalar (integer) cell index
    unsigned int get_scalar_index(const std::vector<unsigned int> 
				  & cell_index) const;

    //! Convert a scalar cell index into a vector cell index
    void get_vector_index(const unsigned int scalar_index,
			  std::vector<unsigned int> & cell_index) const;

  private:
    unsigned int numcells;
    unsigned int spacedims;

    std::vector<unsigned int> blockdims;

    std::vector<T> blockdata;
  };


}  // End namespace Block
}  // End namespace QuickFlash


// Include the function definitions

#include "quickflash_block_blockdata.tcc"


#endif  // QUICKFLASH_BLOCK_BLOCKDATA_HPP
