// C++ program file quickflash_surface_marchingtetrahedrons.cpp

/*
  By Nathan C. Hearn
     July 16, 2008

  Marching Tetrahedrons triangle surface patch generator.

  For a description of the method, see
  http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise
*/


#include "quickflash_surface_marchingtetrahedrons.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include "quickflash_except.hpp"
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
		       Triangle::DoubleVector & interp_coords)
  {
  const unsigned int dims = coords_A.get_dims();

  const double value_A_iso = value_A - isolevel_value;
  const double value_B_iso = value_B - isolevel_value;

  if ((value_A_iso * value_B_iso) > 0.0)
    throw Except("Endpoint values do not enclose the isolevel value", __FILE__,
		 __LINE__);

  // Determine weights for A and B

  const double abs_value_A_B = std::fabs(value_A - value_B);

  double B_weight = 0.5;

  if (abs_value_A_B > 0.0)
    {
    const double abs_value_A_iso = std::fabs(value_A_iso);

    B_weight = std::min(1.0, std::max(0.0, abs_value_A_iso / abs_value_A_B));
    }

  const double A_weight = std::min(1.0, std::max(0.0, 1.0 - B_weight));
  
  // Compose the interpolated vector

  for (unsigned int i = 0; i < dims; i++)
    interp_coords[i] = (A_weight * coords_A[i]) + (B_weight * coords_B[i]);
  }


void get_single_triangle(const Triangle::DoubleVector & apex_coords,
			 const double apex_value,
			 const Triangle::DoubleVector & base_coords_A,
			 const double base_value_A,
			 const Triangle::DoubleVector & base_coords_B,
			 const double base_value_B,
			 const Triangle::DoubleVector & base_coords_C,
			 const double base_value_C,
			 const double isolevel_value, Triangle & triangle)
  {
  Triangle::DoubleVector triangle_node_A;
  Triangle::DoubleVector triangle_node_B;
  Triangle::DoubleVector triangle_node_C;

  find_interp_point(apex_coords, apex_value, base_coords_A, base_value_A,
		    isolevel_value, triangle_node_A);

  find_interp_point(apex_coords, apex_value, base_coords_B, base_value_B,
		    isolevel_value, triangle_node_B);

  find_interp_point(apex_coords, apex_value, base_coords_C, base_value_C,
		    isolevel_value, triangle_node_C);

  triangle.reset(triangle_node_A, triangle_node_B, triangle_node_C);
  }


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
			 const double isolevel_value, Triangle & triangle)
  {
  Triangle::DoubleVector triangle_node_A;
  Triangle::DoubleVector triangle_node_B;
  Triangle::DoubleVector triangle_node_C;

  find_interp_point(apex_coords, apex_value, base_coords_A, base_value_A,
		    isolevel_value, triangle_node_A);

  find_interp_point(apex_coords, apex_value, base_coords_B, base_value_B,
		    isolevel_value, triangle_node_B);

  find_interp_point(base_coords_C_start, base_value_C_start, 
		    base_coords_C_end, base_value_C_end, 
		    isolevel_value, triangle_node_C);

  triangle.reset(triangle_node_A, triangle_node_B, triangle_node_C);
  }


void add_tetra_triangles(const double isolevel_value,
			 const TetrahedronData<const Triangle::DoubleVector *> 
                           & vertex_ptr_list,
                         const TetrahedronData<double> & vertex_values,
			 std::vector<Triangle> & triangle_list)
  {
  // Get the index that describes this tetrahedron

  unsigned int tet_index = 0;

  unsigned int bit_count = 0;

  for (unsigned int i = 0; i < 4; i++)
    if (vertex_values[i] < isolevel_value)
      {
      tet_index += 1 << i;
      bit_count++;
      }

  unsigned int num_triangles = 0;

  if ((tet_index != 0) && (tet_index != 15))
    {
    if ((bit_count % 2) == 1)
      num_triangles = 1;
    else
      num_triangles = 2;
    }

  // Build the triangle list

  triangle_list.resize(num_triangles);

  if (num_triangles > 0)
    {
    // Process the tetrahedron

    if (num_triangles == 1)
      {
      // Node indexes

      unsigned int apex_index = 0;
      unsigned int base_index_A = 0;
      unsigned int base_index_B = 0;
      unsigned int base_index_C = 0;

      switch (tet_index)
	{
	case 1 :
	
	  // Node 0 is below isosurface

	  apex_index = 0;
	
	  // Right-handed normal points away from apex

	  base_index_A = 1;
	  base_index_B = 3;
	  base_index_C = 2;

	  break;

	case 2 :
	
	  // Node 1 is below isosurface

	  apex_index = 1;
	
	  // Right-handed normal points away from apex

	  base_index_A = 0;
	  base_index_B = 2;
	  base_index_C = 3;

	  break;

	case 4 :
	
	  // Node 2 is below isosurface

	  apex_index = 2;
	
	  // Right-handed normal points away from apex

	  base_index_A = 1;
	  base_index_B = 0;
	  base_index_C = 3;

	  break;

	case 8 :
	
	  // Node 3 is below isosurface

	  apex_index = 3;
	
	  // Right-handed normal points away from apex

	  base_index_A = 2;
	  base_index_B = 0;
	  base_index_C = 1;

	  break;

	case 7 :
	
	  // Node 3 is above isosurface
	  
	  apex_index = 3;
	
	  // Right-handed normal points toward apex

	  base_index_A = 0;
	  base_index_B = 2;
	  base_index_C = 1;

	  break;

	case 11 :
	
	  // Node 2 is above isosurface

	  apex_index = 2;
	
	  // Right-handed normal points toward apex

	  base_index_A = 0;
	  base_index_B = 1;
	  base_index_C = 3;

	  break;

	case 13 :
	
	  // Node 1 is above isosurface

	  apex_index = 1;
	
	  // Right-handed normal points toward apex

	  base_index_A = 2;
	  base_index_B = 0;
	  base_index_C = 3;

	  break;

	case 14 :
	
	  // Node 0 is above isosurface

	  apex_index = 0;
	
	  // Right-handed normal points toward apex

	  base_index_A = 3;
	  base_index_B = 1;
	  base_index_C = 2;

	  break;

	default :

	  throw Except("Error in tetrahedron index", __FILE__, __LINE__);
	}

      get_single_triangle(*(vertex_ptr_list[apex_index]), 
			  vertex_values[apex_index],
			  *(vertex_ptr_list[base_index_A]),
			  vertex_values[base_index_A],
			  *(vertex_ptr_list[base_index_B]),
			  vertex_values[base_index_B],
			  *(vertex_ptr_list[base_index_C]),
			  vertex_values[base_index_C],
			  isolevel_value, triangle_list[0]);
      }
    else
      {
      // Node indexes

      unsigned int apex_index_tri1 = 0;
      unsigned int base_index_A_tri1 = 0;
      unsigned int base_index_B_tri1 = 0;
      unsigned int base_index_C_start_tri1 = 0;
      unsigned int base_index_C_end_tri1 = 0;

      unsigned int apex_index_tri2 = 0;
      unsigned int base_index_A_tri2 = 0;
      unsigned int base_index_B_tri2 = 0;
      unsigned int base_index_C_start_tri2 = 0;
      unsigned int base_index_C_end_tri2 = 0;

      switch (tet_index)
	{
	case 3 :
	
	  // Nodes 0 and 1 are below isosurface

	  apex_index_tri1 = 0;
	  apex_index_tri2 = 1;
	
	  // Right-handed normals point away from apexes

	  base_index_A_tri1 = 3;
	  base_index_B_tri1 = 2;
	  base_index_C_start_tri1 = 1;
	  base_index_C_end_tri1 = 2;

	  base_index_A_tri2 = 2;
	  base_index_B_tri2 = 3;
	  base_index_C_start_tri2 = 0;
	  base_index_C_end_tri2 = 3;

	  break;

	case 5 :
	
	  // Nodes 0 and 2 are below isosurface

	  apex_index_tri1 = 0;
	  apex_index_tri2 = 2;
	
	  // Right-handed normals point away from apexes

	  base_index_A_tri1 = 1;
	  base_index_B_tri1 = 3;
	  base_index_C_start_tri1 = 2;
	  base_index_C_end_tri1 = 3;

	  base_index_A_tri2 = 3;
	  base_index_B_tri2 = 1;
	  base_index_C_start_tri2 = 0;
	  base_index_C_end_tri2 = 1;

	  break;

	case 6 :
	
	  // Nodes 1 and 2 are below isosurface

	  apex_index_tri1 = 1;
	  apex_index_tri2 = 2;
	
	  // Right-handed normals point away from apexes

	  base_index_A_tri1 = 3;
	  base_index_B_tri1 = 0;
	  base_index_C_start_tri1 = 2;
	  base_index_C_end_tri1 = 0;

	  base_index_A_tri2 = 0;
	  base_index_B_tri2 = 3;
	  base_index_C_start_tri2 = 1;
	  base_index_C_end_tri2 = 3;

	  break;

	case 9 :
	
	  // Nodes 1 and 2 are above isosurface

	  apex_index_tri1 = 1;
	  apex_index_tri2 = 2;
	
	  // Right-handed normals point toward apexes

	  base_index_A_tri1 = 0;
	  base_index_B_tri1 = 3;
	  base_index_C_start_tri1 = 2;
	  base_index_C_end_tri1 = 0;

	  base_index_A_tri2 = 3;
	  base_index_B_tri2 = 0;
	  base_index_C_start_tri2 = 1;
	  base_index_C_end_tri2 = 3;

	  break;

	case 10 :
	
	  // Nodes 0 and 2 are above isosurface

	  apex_index_tri1 = 0;
	  apex_index_tri2 = 2;
	
	  // Right-handed normals point toward apexes

	  base_index_A_tri1 = 3;
	  base_index_B_tri1 = 1;
	  base_index_C_start_tri1 = 2;
	  base_index_C_end_tri1 = 3;

	  base_index_A_tri2 = 1;
	  base_index_B_tri2 = 3;
	  base_index_C_start_tri2 = 0;
	  base_index_C_end_tri2 = 1;

	  break;

	case 12 :
	
	  // Nodes 0 and 1 are above isosurface

	  apex_index_tri1 = 0;
	  apex_index_tri2 = 1;
	
	  // Right-handed normals point toward apexes

	  base_index_A_tri1 = 2;
	  base_index_B_tri1 = 3;
	  base_index_C_start_tri1 = 1;
	  base_index_C_end_tri1 = 2;

	  base_index_A_tri2 = 3;
	  base_index_B_tri2 = 2;
	  base_index_C_start_tri2 = 0;
	  base_index_C_end_tri2 = 3;

	  break;

	default :

	  throw Except("Error in tetrahedron index", __FILE__, __LINE__);
	}

      // Process triangle 1

      get_single_triangle(*(vertex_ptr_list[apex_index_tri1]), 
			  vertex_values[apex_index_tri1],
			  *(vertex_ptr_list[base_index_A_tri1]),
			  vertex_values[base_index_A_tri1],
			  *(vertex_ptr_list[base_index_B_tri1]),
			  vertex_values[base_index_B_tri1],
			  *(vertex_ptr_list[base_index_C_start_tri1]),
			  vertex_values[base_index_C_start_tri1],
			  *(vertex_ptr_list[base_index_C_end_tri1]),
			  vertex_values[base_index_C_end_tri1],
			  isolevel_value, triangle_list[0]);

      // Process triangle 2

      get_single_triangle(*(vertex_ptr_list[apex_index_tri2]), 
			  vertex_values[apex_index_tri2],
			  *(vertex_ptr_list[base_index_A_tri2]),
			  vertex_values[base_index_A_tri2],
			  *(vertex_ptr_list[base_index_B_tri2]),
			  vertex_values[base_index_B_tri2],
			  *(vertex_ptr_list[base_index_C_start_tri2]),
			  vertex_values[base_index_C_start_tri2],
			  *(vertex_ptr_list[base_index_C_end_tri2]),
			  vertex_values[base_index_C_end_tri2],
			  isolevel_value, triangle_list[1]);
      }
    }
  }


void add_cube_triangles(const double isolevel_value,
			const CubeData<Triangle::DoubleVector>
			  & vertex_positions,
			const CubeData<double> & vertex_values,
			std::list<Triangle> & triangle_list)
  {
  // This is the real interface.  We assume that the list of vertices changes
  // fastest in x and slowest in z.

  // Break up the cube and process the tetrahedra

  TetrahedronData<unsigned int> vertex_index_list;
  vertex_index_list.reset_zero();

  TetrahedronData<const Triangle::DoubleVector *> vertex_ptr_list;
  vertex_ptr_list.reset_zero();

  TetrahedronData<double> vertex_value_list;
  vertex_value_list.reset_zero();

  std::vector<Triangle> triangles;

  // Tetrahedron 1

  vertex_index_list[0] = 0;
  vertex_index_list[1] = 1;
  vertex_index_list[2] = 5;
  vertex_index_list[3] = 2;

  for (unsigned int index = 0; index < 4; index++)
    {
    const unsigned int vertex_index = vertex_index_list[index];

    vertex_ptr_list[index] = &(vertex_positions[vertex_index]);
    vertex_value_list[index] = vertex_values[vertex_index];
    }

  add_tetra_triangles(isolevel_value, vertex_ptr_list, vertex_value_list,
		      triangles);

  for (unsigned int tri_index = 0; tri_index < triangles.size(); tri_index++)
    triangle_list.push_back(triangles[tri_index]);

  // Tetrahedron 2

  vertex_index_list[0] = 0;
  vertex_index_list[1] = 5;
  vertex_index_list[2] = 4;
  vertex_index_list[3] = 2;

  for (unsigned int index = 0; index < 4; index++)
    {
    const unsigned int vertex_index = vertex_index_list[index];

    vertex_ptr_list[index] = &(vertex_positions[vertex_index]);
    vertex_value_list[index] = vertex_values[vertex_index];
    }

  add_tetra_triangles(isolevel_value, vertex_ptr_list, vertex_value_list,
		      triangles);

  for (unsigned int tri_index = 0; tri_index < triangles.size(); tri_index++)
    triangle_list.push_back(triangles[tri_index]);

  // Tetrahedron 3

  vertex_index_list[0] = 4;
  vertex_index_list[1] = 5;
  vertex_index_list[2] = 6;
  vertex_index_list[3] = 2;

  for (unsigned int index = 0; index < 4; index++)
    {
    const unsigned int vertex_index = vertex_index_list[index];

    vertex_ptr_list[index] = &(vertex_positions[vertex_index]);
    vertex_value_list[index] = vertex_values[vertex_index];
    }

  add_tetra_triangles(isolevel_value, vertex_ptr_list, vertex_value_list,
		      triangles);

  for (unsigned int tri_index = 0; tri_index < triangles.size(); tri_index++)
    triangle_list.push_back(triangles[tri_index]);

  // Tetrahedron 4

  vertex_index_list[0] = 2;
  vertex_index_list[1] = 3;
  vertex_index_list[2] = 5;
  vertex_index_list[3] = 7;

  for (unsigned int index = 0; index < 4; index++)
    {
    const unsigned int vertex_index = vertex_index_list[index];

    vertex_ptr_list[index] = &(vertex_positions[vertex_index]);
    vertex_value_list[index] = vertex_values[vertex_index];
    }

  add_tetra_triangles(isolevel_value, vertex_ptr_list, vertex_value_list,
		      triangles);

  for (unsigned int tri_index = 0; tri_index < triangles.size(); tri_index++)
    triangle_list.push_back(triangles[tri_index]);

  // Tetrahedron 5

  vertex_index_list[0] = 6;
  vertex_index_list[1] = 2;
  vertex_index_list[2] = 5;
  vertex_index_list[3] = 7;

  for (unsigned int index = 0; index < 4; index++)
    {
    const unsigned int vertex_index = vertex_index_list[index];

    vertex_ptr_list[index] = &(vertex_positions[vertex_index]);
    vertex_value_list[index] = vertex_values[vertex_index];
    }

  add_tetra_triangles(isolevel_value, vertex_ptr_list, vertex_value_list,
		      triangles);

  for (unsigned int tri_index = 0; tri_index < triangles.size(); tri_index++)
    triangle_list.push_back(triangles[tri_index]);

  // Tetrahedron 6

  vertex_index_list[0] = 1;
  vertex_index_list[1] = 5;
  vertex_index_list[2] = 2;
  vertex_index_list[3] = 3;

  for (unsigned int index = 0; index < 4; index++)
    {
    const unsigned int vertex_index = vertex_index_list[index];

    vertex_ptr_list[index] = &(vertex_positions[vertex_index]);
    vertex_value_list[index] = vertex_values[vertex_index];
    }

  add_tetra_triangles(isolevel_value, vertex_ptr_list, vertex_value_list,
		      triangles);

  for (unsigned int tri_index = 0; tri_index < triangles.size(); tri_index++)
    triangle_list.push_back(triangles[tri_index]);
  }


bool interesting_cube(const double isolevel_value,
		      const CubeData<double> & vertex_values)
  {
  bool ret_val = false;

  const double ref_diff = vertex_values[0] - isolevel_value;

  for (unsigned int corner_index = 1; corner_index < 8; corner_index++)
    {
    const double current_diff = vertex_values[corner_index] - isolevel_value;

    if (!((current_diff * ref_diff) > 0.0))
      ret_val = true;
    }

  return ret_val;
  }


}  // End namespace Surface
}  // End namespace QuickFlash
