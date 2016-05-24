// C++ header file quickflash_geometry_vectorarray.hpp

/*
  By Nathan C. Hearn
     July 16, 2008

  VectorArray container template.
*/


#ifndef QUICKFLASH_GEOMETRY_VECTORARRAY_HPP
#define QUICKFLASH_GEOMETRY_VECTORARRAY_HPP


#include <vector>
#include "quickflash_except.hpp"
#include "quickflash_geometry_vector.hpp"


namespace QuickFlash
{
namespace Geometry
{

template <typename T, unsigned int dims, unsigned int num_vectors>
class VectorArray
  {
  public :

    VectorArray() { reset(); }

    VectorArray(const std::vector< std::vector<T> > & vect_array)
      { reset(vect_array); }

    VectorArray(const std::vector< Vector<T, dims> > & vect_array)
      { reset(vect_array); }

    VectorArray(const VectorArray<T, dims, num_vectors> & source)
      { reset(source); }

    ~VectorArray() { }

    VectorArray<T, dims, num_vectors> 
      & operator=(const VectorArray<T, dims, num_vectors> & source)
      { reset(source); }

    VectorArray<T, dims, num_vectors> 
      & operator=(const std::vector< Vector<T, dims> > & source)
      { reset(source); }

    VectorArray<T, dims, num_vectors> 
      & operator=(const std::vector< std::vector<T> > & source)
      { reset(source); }

    void reset()
      {
      for (unsigned int index = 0; index < num_vectors; index++)
	elems[index].reset();
      }

    void reset(const std::vector< std::vector<T> > & vect_array)
      {
      if (vect_array.size() != num_vectors)
	throw Except("Incompatible number of vectors", __FILE__, __LINE__);

      for (unsigned int index = 0; index < num_vectors; index++)
	elems[index].reset(vect_array[index]);
      }
      
    void reset(const std::vector< Vector<T, dims> > & vect_array)
      {
      if (vect_array.size() != num_vectors)
	throw Except("Incompatible number of vectors", __FILE__, __LINE__);

      for (unsigned int index = 0; index < num_vectors; index++)
	elems[index].reset(vect_array[index]);
      }

    void reset(const VectorArray<T, dims, num_vectors> & source)
      {
      if (&source != this)
	{
	for (unsigned int index = 0; index < num_vectors; index++)
	  elems[index].reset(source[index]);
	}
      }

    unsigned int get_num_vectors() const { return num_vectors; }
    unsigned int get_dims() const { return dims; }

    Vector<T, dims> & operator[](const unsigned int index)
      { return elems[index]; }

    const Vector<T, dims> & operator[](const unsigned int index) const
      { return elems[index]; }

    void copy_vector(const unsigned int index, std::vector<T> & dest)
      { elems[index].copy_elems(dest); }

    void copy_array(std::vector< Vector<T, dims> > & dest)
      {
      dest.resize(num_vectors);

      for (unsigned int index = 0; index < num_vectors; index++)
	dest[index] = elems[index];
      }

    void copy_array(std::vector< std::vector<T> > & dest)
      {
      dest.resize(num_vectors);

      for (unsigned int index = 0; index < num_vectors; index++)
	elems[index].copy_elems(dest[index]);
      }

  private :

    Vector<T, dims> elems[num_vectors];
  };


}  // End namespace Geometry
}  // End namespace QuickFlash


#endif  // QUICKFLASH_GEOMETRY_VECTORARRAY_HPP
