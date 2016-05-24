// C++ program file quickflash_mesh_boundarydefs.cpp

/*
  By Nathan C. Hearn
     August 9, 2007

  Definitions for block boundaries.

  NOTE: There are different gid specifications for Flash2 and Flash 3.

  For Flash2, see: source/database/amr/dBaseDeclarations.F90
  For Flash3, see: source/Simulation/constants.h
*/


#include "quickflash_mesh_boundarydefs.hpp"


namespace QuickFlash
{
namespace Mesh
{

const int NoNeighbor_value_Flash2 = -1;
const int NoNeighbor_value_Flash3 = -1;

const int NoBoundary_value_Flash3 = -10;

const int Boundary_value_min_Flash2 = -50;
const int Boundary_value_max_Flash2 = -20;

const int Boundary_value_min_Flash3 = -50;
const int Boundary_value_max_Flash3 = -20;

const int Reflective_value_Flash2 = -20;
const int Outflow_value_Flash2 = -21;
const int Periodic_value_Flash2 = -22;
const int UserBoundary_value_Flash2 = -23;
const int Hydrostatic_value_Flash2 = -41;

const int ParameshBoundary_value_Flash3 = -20;
const int Reflective_value_Flash3 = -31;
const int Outflow_value_Flash3 = -32;
const int Periodic_value_Flash3 = -35;
const int UserBoundary_value_Flash3 = -38;
const int Isolated_value_Flash3 = -33;
const int Hydrostatic_value_Flash3 = -34;


int get_neighbor_type_flash2(const int gid_value, NeighborType & neighbor_type,
			     BoundaryType & boundary_type)
  {
  int ret_val = -1;  // Error state

  if (gid_value > 0)
    {
    neighbor_type = ValidNeighbor;
    boundary_type = NoBoundary;

    ret_val = 0;  // OK
    }
  else if (gid_value == NoNeighbor_value_Flash2)
    {
    neighbor_type = NoNeighbor;
    boundary_type = NoBoundary;

    ret_val = 0;  // OK
    }
  else if ((gid_value >= Boundary_value_min_Flash2) 
	   && (gid_value <= Boundary_value_max_Flash2))
    {
    neighbor_type = BoundaryNeighbor;

    switch (gid_value)
      {
      case Reflective_value_Flash2 :
	boundary_type = Reflective;
	break;

      case Outflow_value_Flash2 :
	boundary_type = Outflow;
	break;

      case Periodic_value_Flash2 :
	boundary_type = Periodic;
	break;

      case UserBoundary_value_Flash2 :
	boundary_type = UserBoundary;
	break;

      case Hydrostatic_value_Flash2 :
	boundary_type = Hydrostatic;
	break;

      default :
	boundary_type = UnknownBoundary;
      }

    ret_val = 0;  // OK
    }

  return ret_val;
  }


int get_neighbor_type_flash3(const int gid_value, NeighborType & neighbor_type,
			     BoundaryType & boundary_type)
  {
  int ret_val = -1;  // Error state

  if (gid_value > 0)
    {
    neighbor_type = ValidNeighbor;
    boundary_type = NoBoundary;

    ret_val = 0;  // OK
    }
  else if (gid_value == NoNeighbor_value_Flash3)
    {
    neighbor_type = NoNeighbor;
    boundary_type = NoBoundary;

    ret_val = 0;  // OK
    }
  else if (gid_value == NoBoundary_value_Flash3)
    {
    neighbor_type = NoNeighbor;  // REALLY???
    boundary_type = NoBoundary;

    ret_val = 0;  // OK
    }
  else if ((gid_value >= Boundary_value_min_Flash3) 
	   && (gid_value <= Boundary_value_max_Flash3))
    {
    neighbor_type = BoundaryNeighbor;

    switch (gid_value)
      {
      case ParameshBoundary_value_Flash3 :
	boundary_type = UnknownBoundary;
	break;

      case Reflective_value_Flash3 :
	boundary_type = Reflective;
	break;

      case Outflow_value_Flash3 :
	boundary_type = Outflow;
	break;

      case Periodic_value_Flash3 :
	boundary_type = Periodic;
	break;

      case UserBoundary_value_Flash3 :
	boundary_type = UserBoundary;
	break;

      case Isolated_value_Flash3 :
	boundary_type = Isolated;
	break;

      case Hydrostatic_value_Flash3 :
	boundary_type = Hydrostatic;
	break;

      default :
	boundary_type = UnknownBoundary;
      }

    ret_val = 0;  // OK
    }

  return ret_val;
  }


}  // End namespace Mesh
}  // End namespace QuickFlash
