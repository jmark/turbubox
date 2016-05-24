// C++ header file quickflash_surface_triangle.hpp

/*
  By Nathan C. Hearn
     July 16, 2008

  Triangle patch description class.

  NOTE: Outward normal vectors are defined using right-handed ordering for
  vertices.

  NOTE: For now, this class is designed for 3D operation only.
*/


#ifndef QUICKFLASH_SURFACE_TRIANGLE_HPP
#define QUICKFLASH_SURFACE_TRIANGLE_HPP


#include <vector>
#include "quickflash_geometry_vector.hpp"


namespace QuickFlash
{
namespace Surface
{

// Class Triangle

class Triangle
  {
  public :

    typedef Geometry::Vector<double, 3> DoubleVector;

  public :

    Triangle() : normal_vect() { reset(); }

    Triangle(const std::vector<double> & vertex_1, 
	     const std::vector<double> & vertex_2,
	     const std::vector<double> & vertex_3) :
      normal_vect()
      { reset(vertex_1, vertex_2, vertex_3); }

    Triangle(const DoubleVector & vertex_1, 
	     const DoubleVector & vertex_2,
	     const DoubleVector & vertex_3) :
      normal_vect()
      { reset(vertex_1, vertex_2, vertex_3); }

    Triangle(const std::vector< std::vector<double> > & vertex_list) :
      normal_vect()
      { reset(vertex_list); }

    Triangle(const std::vector<DoubleVector> & vertex_list) :
      normal_vect()
      { reset(vertex_list); }

    Triangle(const Triangle & source) :
      normal_vect()
      { reset(source); }

    ~Triangle() { }

    Triangle & operator=(const Triangle & source)
      { 
      reset(source); 
      return *this;
      }

    void reset();

    void reset(const std::vector<double> & vertex_1, 
	       const std::vector<double> & vertex_2,
	       const std::vector<double> & vertex_3);

    void reset(const DoubleVector & vertex_1, 
	       const DoubleVector & vertex_2,
	       const DoubleVector & vertex_3);

    void reset(const std::vector< std::vector<double> > & vertex_list);

    void reset(const std::vector<DoubleVector> & vertex_list);

    void reset(const Triangle & source);

    DoubleVector & get_vertex(const unsigned int index)
      { return vertex_coords[index]; }

    const DoubleVector & get_vertex(const unsigned int index) const
      { return vertex_coords[index]; }

    DoubleVector & operator[](const unsigned int index)
      { return get_vertex(index); }

    const DoubleVector & operator[](const unsigned int index) const
      { return get_vertex(index); }

    void copy_vertex(const unsigned int index, std::vector<double> & dest) 
      const
      { vertex_coords[index].copy_elems(dest); }

    double get_surface_area() const;

    const DoubleVector & get_outward_unit_normal() const
      {
      // DANGER: NOT SAFE FROM CHANGES TO ELEMENTS!!!!!

      if (!(normal_vect_computed))
	compute_outward_normal_vect();

      return normal_vect;
      }

    void copy_unit_normal(std::vector<double> & dest) const
      { get_outward_unit_normal().copy_elems(dest); }


  private :

    void compute_outward_normal_vect() const;


  private :
    
    static const unsigned int dims;

    DoubleVector vertex_coords[3];

    mutable bool normal_vect_computed;

    mutable DoubleVector normal_vect;
  };


}  // End namespace Surface
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SURFACE_TRIANGLE_HPP
