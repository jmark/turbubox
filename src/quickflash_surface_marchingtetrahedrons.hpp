// C++ header file quickflash_surface_marchingtetrahedrons.hpp

/*
  By Nathan C. Hearn
     July 16, 2008

  Marching Tetrahedrons triangle surface patch generator.

  For a description of the method, see
  http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise
*/


#ifndef QUICKFLASH_SURFACE_MARCHINGTETRAHEDRONS_HPP
#define QUICKFLASH_SURFACE_MARCHINGTETRAHEDRONS_HPP


#include <vector>
#include <list>
#include "quickflash_surface_triangle.hpp"
#include "quickflash_surface_cube.hpp"
#include "quickflash_surface_tetrahedron.hpp"


namespace QuickFlash
{
namespace Surface
{

void find_interp_point(const Triangle::DoubleVector & coords_A,
		       const double value_A,
		       const Triangle::DoubleVector & coords_B,
		       const double value_B,
		       const double isolevel_value,
		       Triangle::DoubleVector & interp_coords);


void get_single_triangle(const Triangle::DoubleVector & apex_coords,
			 const double apex_value,
			 const Triangle::DoubleVector & base_coords_A,
			 const double base_value_A,
			 const Triangle::DoubleVector & base_coords_B,
			 const double base_value_B,
			 const Triangle::DoubleVector & base_coords_C,
			 const double base_value_C,
			 const double isolevel_value, 
			 Triangle & triangle);


void get_single_triangle(const Triangle::DoubleVector & apex_coords,
			 const double apex_value,
			 const Triangle::DoubleVector & base_coords_A,
			 const double base_value_A,
			 const Triangle::DoubleVector & base_coords_B,
			 const double base_value_B,
			 const Triangle::DoubleVector & base_coords_C_start,
			 const double base_value_C_start,
			 const Triangle::DoubleVector & base_coords_C_end,
			 const double base_value_C_end,
			 const double isolevel_value, 
			 Triangle & triangle);


void add_tetra_triangles(const double isolevel_value,
			 const TetrahedronData<const Triangle::DoubleVector *> 
			   & vertex_ptr_list,
			 const TetrahedronData<double> & vertex_values,
			 std::vector<Triangle> & triangle_list);


void add_cube_triangles(const double isolevel_value,
			const CubeData<Triangle::DoubleVector> 
			  & vertex_positions,
			const CubeData<double> & vertex_values,
			std::list<Triangle> & triangle_list);


bool interesting_cube(const double isolevel_value,
		      const CubeData<double> & vertex_values);


}  // End namespace Surface
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SURFACE_MARCHINGTETRAHEDRONS_HPP
