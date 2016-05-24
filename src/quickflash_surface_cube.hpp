// C++ header file quickflash_surface_cube.hpp

/*
    By Nathan C. Hearn
     May 8, 2008

  Container for cube vertex data.

  NOTE: For 3D index vectors, data is stored with the first index changing 
        the fastest, and the last index changing the slowest.
*/


#ifndef QUICKFLASH_SURFACE_CUBE_HPP
#define QUICKFLASH_SURFACE_CUBE_HPP


namespace QuickFlash
{
namespace Surface
{

// Class CubeData

template <typename T>
class CubeData
  {
  public :

    CubeData() { }

    CubeData(const CubeData<T> & source) { reset(source); }

    ~CubeData() { }

    CubeData<T> & operator=(const CubeData<T> & source)
      {
      reset(source);
      return *this;
      }

    void reset()
      {
      for (unsigned int corner_index = 0; corner_index < 8; corner_index++)
	elements[corner_index].reset();
      }

    void reset(const T & value)
      {
      for (unsigned int corner_index = 0; corner_index < 8; corner_index++)
	elements[corner_index] = value;
      }

    void reset(const CubeData<T> & source)
      {
      if (&source != this)
	{
	for (unsigned int corner_index = 0; corner_index < 8; corner_index++)
	  elements[corner_index] = source.elements[corner_index];
	}
      }

    void reset_zero()
      {
      for (unsigned int corner_index = 0; corner_index < 8; corner_index++)
	elements[corner_index] = static_cast<T>(0);
      }

    T & operator[](const unsigned int index)
      { return get_elem(index); }

    const T & operator[](const unsigned int index) const
      { return get_elem(index); }

    T & get_elem(const unsigned int index) { return elements[index]; }

    const T & get_elem(const unsigned int index) const 
      { return elements[index]; }
			 
    T & get_elem(const unsigned int i, const unsigned int j, 
		 const unsigned int k)
      { return get_elem(i + (2 * (j + (2 * k)))); }

    const T & get_elem(const unsigned int i, const unsigned int j, 
		       const unsigned int k) const
      { return get_elem(i + (2 * (j + (2 * k)))); }

  private :

    T elements[8];
  };


}  // End namespace Surface
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SURFACE_CUBE_HPP
