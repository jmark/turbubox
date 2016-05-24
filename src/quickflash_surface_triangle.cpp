// C++ program file quickflash_surface_triangle.cpp

/*
  By Nathan C. Hearn
     July 16, 2008

  Triangle patch description class.

  NOTE: For now, this class is designed for 3D operation only.
*/


#include "quickflash_surface_triangle.hpp"
#include <vector>
#include <cmath>
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Surface
{

// Class Triangle

const unsigned int Triangle::dims = 3;


void Triangle::reset()
  {
  for (unsigned int vert_index = 0; vert_index < 3; vert_index++)
    vertex_coords[vert_index].reset();

  normal_vect_computed = false;

  normal_vect.reset();
  }


void Triangle::reset(const std::vector<double> & vertex_1, 
		     const std::vector<double> & vertex_2,
		     const std::vector<double> & vertex_3)
  {
  // Store the data

  vertex_coords[0] = vertex_1;
  vertex_coords[1] = vertex_2;
  vertex_coords[2] = vertex_3;

  normal_vect_computed = false;
  }


void Triangle::reset(const Triangle::DoubleVector & vertex_1,
		     const Triangle::DoubleVector & vertex_2,
		     const Triangle::DoubleVector & vertex_3)
  {
  // Store the data

  vertex_coords[0] = vertex_1;
  vertex_coords[1] = vertex_2;
  vertex_coords[2] = vertex_3;

  normal_vect_computed = false;
  }
  

void Triangle::reset(const std::vector< std::vector<double> > & vertex_list)
  {
  // Check the sizes

  if (vertex_list.size() != 3)
    throw Except("Three vertices required", __FILE__, __LINE__);

  // Store the data

  for (unsigned int vert_index = 0; vert_index < 3; vert_index++)
    vertex_coords[vert_index] = vertex_list[vert_index];

  normal_vect_computed = false;
  }


void Triangle::reset(const std::vector<Triangle::DoubleVector> & vertex_list)
  {
  // Check the sizes

  if (vertex_list.size() != 3)
    throw Except("Three vertices required", __FILE__, __LINE__);

  // Store the data

  for (unsigned int vert_index = 0; vert_index < 3; vert_index++)
    vertex_coords[vert_index] = vertex_list[vert_index];

  normal_vect_computed = false;
  }


void Triangle::reset(const Triangle & source)
  {
  if (&source != this)
    {
    for (unsigned int vert_index = 0; vert_index < 3; vert_index++)
      vertex_coords[vert_index] = source.vertex_coords[vert_index];

    normal_vect_computed = source.normal_vect_computed;

    normal_vect = source.normal_vect;
    }
  }


double Triangle::get_surface_area() const
  {
  // REALLY ... THIS SHOULD BE CHANGED TO USE CROSS-PRODUCTS!

  double ret_val = 0.0;  // Area is zero unless found otherwise

  const DoubleVector & origin = vertex_coords[0];

  // Let v0--v1 be the base of the triangle

  DoubleVector base = vertex_coords[1];

  base -= origin;

  const double base_sq_len = base.get_sq_length();

  if (base_sq_len > 0.0)
    {
    // Use the angle between the base and v0--v2 (adj_side) to get the height

    DoubleVector adj_side = vertex_coords[2];

    adj_side -= origin;

    const double adj_sq_len = adj_side.get_sq_length();

    if (adj_sq_len > 0.0)
      {
      // Use the dot product to get the cosine of the angle

      const double dot_prod = adj_side.dot_prod(base);

      const double adj_base_sq_len_prod = adj_sq_len * base_sq_len;

      // dot_prod^2 = cos(angle)^2 * base_len^2 * adj_len^2

      const double sq_cos_angle = (dot_prod * dot_prod) / adj_base_sq_len_prod;

      // Area is base * height / 2
      // ... height is |v0--v2| * sin(angle)

      const double sq_sin_angle = 1.0 - sq_cos_angle;

      if (sq_sin_angle > 0.0)
	ret_val = 0.5 * std::sqrt(adj_base_sq_len_prod * sq_sin_angle);
      }
    }

  return ret_val;
  }


void Triangle::compute_outward_normal_vect() const
  {
  // Assume a right-handed order of vertexes 
  //   (i.e., normal = (v1 - v0) x (v2 - v1)

  const Triangle::DoubleVector & v0 = vertex_coords[0];
  const Triangle::DoubleVector & v1 = vertex_coords[1];
  const Triangle::DoubleVector & v2 = vertex_coords[2];

  Triangle::DoubleVector v1_v0;

  for (unsigned int i = 0; i < dims; i++)
    v1_v0[i] = v1[i] - v0[i];

  Triangle::DoubleVector v2_v1;

  for (unsigned int i = 0; i < dims; i++)
    v2_v1[i] = v2[i] - v1[i];

  normal_vect[0] = (v1_v0[1] * v2_v1[2]) - (v1_v0[2] * v2_v1[1]);
  normal_vect[1] = (v1_v0[2] * v2_v1[0]) - (v1_v0[0] * v2_v1[2]);
  normal_vect[2] = (v1_v0[0] * v2_v1[1]) - (v1_v0[1] * v2_v1[0]);

  // Get the unit normal

  double vect_sq_len = 0.0;

  for (unsigned int i = 0; i < dims; i++)
    {
    const double elem = normal_vect[i];

    vect_sq_len += elem * elem;
    }

  if (!(vect_sq_len > 0.0))
    throw Except("Error computing unit normal", __FILE__, __LINE__);

  const double vect_len = std::sqrt(vect_sq_len);

  const double oneOver_vect_len = 1.0 / vect_len;

  for (unsigned int i = 0; i < dims; i++)
    normal_vect[i] *= oneOver_vect_len;

  normal_vect_computed = true;
  }


}  // End namespace Surface
}  // End namespace QuickFlash
