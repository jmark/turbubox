// C++ header file quickflash_geometry_vector.hpp

/*
  By Nathan C. Hearn
     July 16, 2008

  Vector container template.
*/


#ifndef QUICKFLASH_GEOMETRY_VECTOR_HPP
#define QUICKFLASH_GEOMETRY_VECTOR_HPP


#include <vector>
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Geometry
{

template <typename T, unsigned int dims>
class Vector
  {
  public :

    Vector() { reset(); }

    Vector(const std::vector<T> & vect_elems) { reset(vect_elems); }

    Vector(const Vector<T, dims> & source) { reset(source); }

    ~Vector() { }

    Vector<T, dims> & operator=(const Vector<T, dims> & source)
      {
      reset(source);
      return *this;
      }

    Vector<T, dims> & operator=(const std::vector<T> & source)
      {
      reset(source);
      return *this;
      }

    Vector<T, dims> & operator+=(const Vector<T, dims> & operand)
      {
      for (unsigned int i = 0; i < dims; i++)
	elems[i] += operand.elems[i];

      return *this;
      }

    Vector<T, dims> & operator-=(const Vector<T, dims> & operand)
      {
      for (unsigned int i = 0; i < dims; i++)
	elems[i] -= operand.elems[i];

      return *this;
      }

    void reset()
      {
      for (unsigned int i = 0; i < dims; i++)
	elems[i] = static_cast<T>(0);
      }

    void reset(const std::vector<T> & vect_elems)
      {
      if (vect_elems.size() != dims)
	throw Except("Incompatible vector length", __FILE__, __LINE__);

      for (unsigned int i = 0; i < dims; i++)
	elems[i] = vect_elems[i];
      }

    void reset(const Vector<T, dims> & source)
      {
      if (&source != this)
	{
	for (unsigned int i = 0; i < dims; i++)
	  elems[i] = source.elems[i];
	}
      }

    unsigned int get_dims() const { return dims; }

    T & operator[](const unsigned int index)
      { return elems[index]; }

    const T & operator[](const unsigned int index) const
      { return elems[index]; }

    void copy_elems(std::vector<T> & dest) const
      {
      dest.resize(dims);

      for (unsigned int i = 0; i < dims; i++)
	dest[i] = elems[i];
      }

    double get_sq_length() const
      {
      double sq_len = 0.0;

      for (unsigned int i = 0; i < dims; i++)
	{
	const double current_elem = elems[i];

	sq_len += current_elem * current_elem;
	}

      return sq_len;
      }

    double dot_prod(const Vector<T, dims> & operand)
      {
      double product = 0.0;

      for (unsigned int i = 0; i < dims; i++)
	product += elems[i] * operand.elems[i];

      return product;
      }

  private :

    T elems[dims];
  };


}  // End namespace Geometry
}  // End namespace QuickFlash


#endif  // QUICKFLASH_GEOMETRY_VECTOR_HPP
