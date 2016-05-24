// C++ header file quickflash_geometry.hpp

/*
  By Nathan C. Hearn
     June 1, 2007

  Geometry definitions and utilities for QuickFlash.
*/


#ifndef QUICKFLASH_GEOMETRY_HPP
#define QUICKFLASH_GEOMETRY_HPP


#include <string>
#include <vector>


namespace QuickFlash
{
namespace Geometry
{

// Basic definitions

extern const double two_pi;
extern const double oneOver_two_pi;
extern const double fourThirds_pi;


// Mesh geometry types

enum GeometryType { Cartesian, Polar, Cylindrical, Spherical};

extern const char Cartesian_Label[];
extern const char Polar_Label[];
extern const char Cylindrical_Label[];
extern const char Spherical_Label[];

int get_geometry_type(const std::string & label, GeometryType & geometry_type);

const char * get_geometry_label(const GeometryType geometry_type);
void get_geometry_label(const GeometryType geometry_type, std::string & label);


// Flash coordinate indexes

extern const unsigned int Polar_Radius_Index;
extern const unsigned int Polar_Azimuth_Index;

extern const unsigned int Cylindrical_Radius_Index;
extern const unsigned int Cylindrical_zAxis_Index;
extern const unsigned int Cylindrical_Azimuth_Index;

extern const unsigned int Spherical_Radius_Index;
extern const unsigned int Spherical_Polar_Index;
extern const unsigned int Spherical_Azimuth_Index;


// Physical volume of logical cubic domain

int get_domain_volume_phys(const GeometryType geometry_type,
			   const std::vector<double> & min_domain_bounds,
			   const std::vector<double> & max_domain_bounds,
			   double & domain_volume);

int get_position_cartesian(const GeometryType geometry_type,
			   const std::vector<double> & logical_coords,
			   std::vector<double> & cartesian_coords);

int get_origin_sqDist(const GeometryType geometry_type,
		    const std::vector<double> & logical_coords,
		    double & origin_sqDistance);

int get_polar_axis_sqDist(const GeometryType geometry_type,
			  const std::vector<double> & logical_coords,
			  double & polar_axis_sqDistance);


}  // End namespace Geometry
}  // End namespace QuickFlash


#endif  // QUICKFLASH_GEOMETRY_HPP
