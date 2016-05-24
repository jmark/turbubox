// C++ header file quickflash_mesh_boundarydefs.hpp

/*
  By Nathan C. Hearn
     August 9, 2007

  Definitions for block neighbors.
*/


#ifndef QUICKFLASH_MESH_BOUNDARYDEFS_HPP
#define QUICKFLASH_MESH_BOUNDARYDEFS_HPP


namespace QuickFlash
{
namespace Mesh
{

enum NeighborType { NoNeighbor, BoundaryNeighbor, ValidNeighbor };
enum BoundaryType { NoBoundary, UnknownBoundary, Reflective, 
		    Outflow, Periodic, Hydrostatic, Isolated, UserBoundary };

extern const int NoNeighbor_value_Flash2;  // No neighbor at same refinement
extern const int NoNeighbor_value_Flash3;

extern const int NoBoundary_value_Flash3;

extern const int Boundary_value_min_Flash2;
extern const int Boundary_value_max_Flash2;

extern const int Boundary_value_min_Flash3;
extern const int Boundary_value_max_Flash3;

extern const int Reflective_value_Flash2;
extern const int Outflow_value_Flash2;
extern const int Periodic_value_Flash2;
extern const int UserBoundary_value_Flash2;
extern const int Hydrostatic_value_Flash2;

extern const int ParameshBoundary_value_Flash3;
extern const int Reflective_value_Flash3;
extern const int Outflow_value_Flash3;
extern const int Periodic_value_Flash3;
extern const int UserBoundary_value_Flash3;
extern const int Isolated_value_Flash3;
extern const int Hydrostatic_value_Flash3;


int get_neighbor_type_flash2(const int gid_value, 
			     NeighborType & neighbor_type, 
			     BoundaryType & boundary_type);

int get_neighbor_type_flash3(const int gid_value, 
			     NeighborType & neighbor_type, 
			     BoundaryType & boundary_type);


}  // End namespace Mesh
}  // End namespace QuickFlash


#endif  // QUICKFLASH_MESH_BOUNDARYDEFS_HPP
