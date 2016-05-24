// C++ auxiliary template file quickflash_block_blockdata.tcc

/*
  By Nathan C. Hearn
     May 23, 2007

  For inclusion only by quickflash_block_blockdata.hpp.
*/


#ifndef QUICKFLASH_BLOCK_BLOCKDATA_TCC
#define QUICKFLASH_BLOCK_BLOCKDATA_TCC


#include "quickflash_block_blockdata.hpp"


namespace QuickFlash
{
namespace Block
{

// Template Class BlockData


template<typename T>
void BlockData<T>::get_minmax_values(T & min_value, T & max_value)
  const
  {
  // This function requires the less-than operator (<) for type T

  if (numcells < 1)
    throw Except("Block uninitialized", __FILE__, __LINE__);

  T min_val = blockdata[0];
  T max_val = min_val;

  for (unsigned int index = 1; index < numcells; index++)
    {
    const T & value = blockdata[index];

    if (value < min_val)
      min_val = value;
    else if (!(value < max_val))
      max_val = value;
    }

  min_value = min_val;
  max_value = max_val;
  }


template<typename T>
void BlockData<T>::reset(const unsigned int space_dims, 
			 const unsigned int block_dims[])
  {
  if (space_dims < 1)
    throw Except("No spatial dimensions", __FILE__, __LINE__);

  spacedims = space_dims;

  blockdims.resize(spacedims);
  numcells = 1;

  for (unsigned int i = 0; i < space_dims; i++)
    {
    const unsigned int axis_dims = block_dims[i];

    if (axis_dims < 1)
      throw Except("Improper number of cells", __FILE__, __LINE__);

    blockdims[i] = axis_dims;

    numcells *= axis_dims;
    }

  blockdata.resize(numcells);  // Does not clear data!!!
  }


template<typename T>
void BlockData<T>::reset(const std::vector<unsigned int> & block_dims)
  {
  spacedims = block_dims.size();

  if (spacedims < 1)
    throw Except("No spatial dimensions", __FILE__, __LINE__);

  blockdims = block_dims;
  numcells = 1;

  for (unsigned int i = 0; i < spacedims; i++)
    {
    const unsigned int axis_dims = block_dims[i];

    if (axis_dims < 1)
      throw Except("Improper number of cells", __FILE__, __LINE__);

    numcells *= axis_dims;
    }

  blockdata.resize(numcells);  // Does not clear data!!!
  }


template<typename T>
void BlockData<T>::set_data_reverse(const unsigned int num_elems, 
				    const T block_data[])
  {
  if (num_elems != blockdata.size())
    throw Except("Incompatible data vector", __FILE__, __LINE__);

  // Correct for Fortran-style data ordering

  VectorCounter counter(blockdims);

  while (counter.in_bounds())
    {
    const unsigned int fortran_index = counter.get_scalar_index_reverse();
    const unsigned int local_index = counter.get_scalar_index();

    blockdata[local_index] = block_data[fortran_index];

    counter.increment();
    }
  }


template<typename T>
unsigned int BlockData<T>::get_scalar_index(const std::vector<unsigned int> 
					    & cell_index) const
  {
  unsigned int scalar_index = cell_index[0];

  for (unsigned int i = 1; i < spacedims; i++)
    {
    scalar_index *= blockdims[i];
    scalar_index += cell_index[i];
    }

  return scalar_index;
  }


template<typename T>
void BlockData<T>::get_vector_index(const unsigned int scalar_index,
				    std::vector<unsigned int> & cell_index) 
  const
  {
  unsigned int temp_scalar_index = scalar_index;

  cell_index.resize(spacedims);

  for (unsigned int i = (spacedims - 1); i > 0; i--)
    {
    const unsigned int current_axis_dim = blockdims[i];

    const unsigned int current_axis_size
      = temp_scalar_index % current_axis_dim;

    cell_index[i] = current_axis_size;

    temp_scalar_index /= current_axis_dim;
    }

  cell_index[0] = temp_scalar_index;
  }


}  // End namespace Block
}  // End namespace QuickFlash


#endif  // QUICKFLASH_BLOCK_BLOCKDATA_TCC
