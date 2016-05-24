// C++ header file quickflash_surface_tetrahedron.hpp

/*
    By Nathan C. Hearn
     July 16, 2008

  Container for tetrahedron vertex data.

  NOTE: For 3D index vectors, data is stored with the first index changing 
        the fastest, and the last index changing the slowest.
*/


#ifndef QUICKFLASH_SURFACE_TETRAHEDRON_HPP
#define QUICKFLASH_SURFACE_TETRAHEDRON_HPP


namespace QuickFlash
{
namespace Surface
{

// Class TetrahedronData

template <typename T>
class TetrahedronData
  {
  public :

    TetrahedronData() { }

    TetrahedronData(const TetrahedronData<T> & source) 
      { reset(source); }

    ~TetrahedronData() { }

    TetrahedronData<T> & operator=(const TetrahedronData<T> & source)
      {
      reset(source);
      return *this;
      }

    void reset()
      {
      for (unsigned int corner_index = 0; corner_index < 4; corner_index++)
	elements[corner_index].reset();
      }

    void reset(const T & value)
      {
      for (unsigned int corner_index = 0; corner_index < 4; corner_index++)
	elements[corner_index] = value;
      }

    void reset(const TetrahedronData<T> & source)
      {
      if (&source != this)
	{
	for (unsigned int corner_index = 0; corner_index < 4; corner_index++)
	  elements[corner_index] = source.elements[corner_index];
	}
      }

    void reset_zero()
      {
      for (unsigned int corner_index = 0; corner_index < 4; corner_index++)
	elements[corner_index] = static_cast<T>(0);
      }

    T & operator[](const unsigned int index)
      { return get_elem(index); }

    const T & operator[](const unsigned int index) const
      { return get_elem(index); }

    T & get_elem(const unsigned int index) { return elements[index]; }

    const T & get_elem(const unsigned int index) const 
      { return elements[index]; }

  private :

    T elements[4];
  };


}  // End namespace Surface
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SURFACE_TETRAHEDRON_HPP
