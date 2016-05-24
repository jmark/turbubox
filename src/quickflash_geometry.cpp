// C++ program file quickflash_geometry.cpp

/*
  By Nathan C. Hearn
     June 1, 2007

  Geometry definitions and utilities for QuickFlash.
*/


#include "quickflash_geometry.hpp"
#include <string>
#include <vector>
#include <math.h>


namespace QuickFlash
{
namespace Geometry
{

// Basic definitions

const double two_pi = 2.0 * (M_PI);
const double oneOver_two_pi = 1.0 / two_pi;
const double fourThirds_pi = (4.0 / 3.0) * (M_PI);


// Mesh geometry types

const char Cartesian_Label[] = "cartesian";
const char Polar_Label[] = "polar";
const char Cylindrical_Label[] = "cylindrical";
const char Spherical_Label[] = "spherical";


int get_geometry_type(const std::string & label, GeometryType & geometry_type)
  {
  int ret_val = -1;  // Error state

  if (label == Cartesian_Label)
    {
    geometry_type = Cartesian;
    ret_val = 0;  // OK
    }
  else if (label == Polar_Label)
    {
    geometry_type = Polar;
    ret_val = 0;  // OK
    }
  else if (label == Cylindrical_Label)
    {
    geometry_type = Cylindrical;
    ret_val = 0;  // OK
    }
  else if (label == Spherical_Label)
    {
    geometry_type = Spherical;
    ret_val = 0;  // OK
    }

  return ret_val;
  }


const char * get_geometry_label(const GeometryType geometry_type)
  {
  const char * label_ptr = 0;

  switch (geometry_type)
    {
    case Cartesian :

      label_ptr = Cartesian_Label;
      break;

    case Polar :

      label_ptr = Polar_Label;
      break;

    case Cylindrical :

      label_ptr = Cylindrical_Label;
      break;

    case Spherical :

      label_ptr = Spherical_Label;
      break;
    }

  return label_ptr;
  }


void get_geometry_label(const GeometryType geometry_type, std::string & label)
  { label = get_geometry_label(geometry_type); }


// Flash coordinate indexes

const unsigned int Polar_Radius_Index = 0;
const unsigned int Polar_Azimuth_Index = 1;

const unsigned int Cylindrical_Radius_Index = 0;
const unsigned int Cylindrical_zAxis_Index = 1;
const unsigned int Cylindrical_Azimuth_Index = 2;

const unsigned int Spherical_Radius_Index = 0;
const unsigned int Spherical_Polar_Index = 1;
const unsigned int Spherical_Azimuth_Index = 2;


// Physical volume of logical cubic domain

int get_domain_volume_phys(const GeometryType geometry_type,
			   const std::vector<double> & min_logical_bounds,
			   const std::vector<double> & max_logical_bounds,
			   double & physical_volume)
  {
  int ret_val = -1;  // Error state

  // ERROR CHECKING???

  const unsigned int dims = min_logical_bounds.size();

  if ((dims >= 1) && (dims <= 3))
    {
    double volume = 0.0;

    switch (geometry_type)
      {
      case Cartesian :
	{
	volume = 1.0;

	for (unsigned int i = 0; i < dims; i++)
	  volume *= (max_logical_bounds[i] - min_logical_bounds[i]);

	ret_val = 0;  // OK
	}
	break;

      case Polar :
	{
	// Only valid for 1 or 2 dimensions

	const double min_r = min_logical_bounds[Polar_Radius_Index];
	const double max_r = max_logical_bounds[Polar_Radius_Index];

	volume = (M_PI) * ((max_r * max_r) - (min_r * min_r));

	if (dims > 1)
	  {
	  const double min_azim = min_logical_bounds[Polar_Azimuth_Index];
	  const double max_azim = max_logical_bounds[Polar_Azimuth_Index];

	  volume *= oneOver_two_pi * (max_azim - min_azim);
	  }

	ret_val = 0;  // OK
	}
	break;

      case Cylindrical :
	{
	const double min_r = min_logical_bounds[Cylindrical_Radius_Index];
	const double max_r = max_logical_bounds[Cylindrical_Radius_Index];

	volume = (M_PI) * ((max_r * max_r) - (min_r * min_r));

	if (dims > 1)
	  {
	  const double min_z = min_logical_bounds[Cylindrical_zAxis_Index];
	  const double max_z = max_logical_bounds[Cylindrical_zAxis_Index];

	  volume *= (max_z - min_z);
	  }

	if (dims > 2)
	  {
	  const double min_azim 
	    = min_logical_bounds[Cylindrical_Azimuth_Index];
	  const double max_azim 
	    = max_logical_bounds[Cylindrical_Azimuth_Index];

	  volume *= oneOver_two_pi * (max_azim - min_azim);
	  }

	ret_val = 0;  // OK
	}
	break;

      case Spherical :
	{
	const double min_r = min_logical_bounds[Spherical_Radius_Index];
	const double max_r = max_logical_bounds[Spherical_Radius_Index];

	volume = fourThirds_pi 
	  * ((max_r * max_r * max_r) - (min_r * min_r * min_r));

	if (dims > 1)
	  {
	  const double min_polar = min_logical_bounds[Spherical_Polar_Index];
	  const double max_polar = max_logical_bounds[Spherical_Polar_Index];

	  volume *= 0.5 * (cos(min_polar) - cos(max_polar));
	  }

	if (dims > 2)
	  {
	  const double min_azim = min_logical_bounds[Spherical_Azimuth_Index];
	  const double max_azim = max_logical_bounds[Spherical_Azimuth_Index];

	  volume *= oneOver_two_pi * (max_azim - min_azim);
	  }

	ret_val = 0;  // OK
	}
	break;
      }

    physical_volume = volume;
    }

  return ret_val;
  }


int get_position_cartesian(const GeometryType geometry_type,
			    const std::vector<double> & logical_coords,
			    std::vector<double> & cartesian_coords)
  {
  int ret_val = -1;  // Error state

  const unsigned int dims = logical_coords.size();

  if ((dims >= 1) && (dims <= 3))
    {
    cartesian_coords.resize(dims);

    for (unsigned int i = 0; i < dims; i++)
      cartesian_coords[i] = 0.0;

    switch (geometry_type)
      {
      case Cartesian :
	{
	cartesian_coords = logical_coords;
	ret_val = 0;  // OK
	}
	break;

      case Polar :
	{
	// Only valid for one or two dimensions

	const double r = logical_coords[Polar_Radius_Index];

	if (dims < 2)
	  cartesian_coords[0] = r;
	else
	  {
	  const double azim = logical_coords[Polar_Azimuth_Index];

	  cartesian_coords[0] = r * cos(azim);
	  cartesian_coords[1] = r * sin(azim);
	  }

	ret_val = 0;  // OK
	}
	break;

      case Cylindrical :
	{
	const double xy_dist = logical_coords[Cylindrical_Radius_Index];

	if (dims < 2)
	  cartesian_coords[0] = xy_dist;
	else
	  {
	  const double z = logical_coords[Cylindrical_zAxis_Index];

	  if (dims < 3)
	    {
	    cartesian_coords[0] = xy_dist;
	    cartesian_coords[1] = z;  // Just treat z as y in 2D
	    }
	  else
	    {
	    const double azim = logical_coords[Cylindrical_Azimuth_Index];

	    cartesian_coords[0] = xy_dist * cos(azim);
	    cartesian_coords[1] = xy_dist * sin(azim);
	    cartesian_coords[2] = z;
	    }
	  }

	ret_val = 0;  // OK
	}
	break;

      case Spherical :
	{
	const double r = logical_coords[Spherical_Radius_Index];

	if (dims < 2)
	  cartesian_coords[0] = r;
	else
	  {
	  const double polar = logical_coords[Spherical_Polar_Index];

	  const double xy_dist = r * sin(polar);
	  const double z = r * cos(polar);

	  if (dims < 3)
	    {
	    cartesian_coords[0] = xy_dist;  // Treat as point in x-y plane
	    cartesian_coords[1] = z;        // Treat z as y in 2D
	    }
	  else
	    {
	    const double azim = logical_coords[Spherical_Azimuth_Index];

	    cartesian_coords[0] = xy_dist * cos(azim);
	    cartesian_coords[1] = xy_dist * sin(azim);
	    cartesian_coords[2] = z;
	    }
	  }

	ret_val = 0;  // OK
	}
	break;
      }
    }

  return ret_val;
  }


int get_origin_sqDist(const GeometryType geometry_type,
		    const std::vector<double> & logical_coords,
		    double & origin_sqDistance)
  {
  int ret_val = -1;  // Error state

  const unsigned int dims = logical_coords.size();

  if ((dims >= 1) && (dims <= 3))
    {
    double sq_dist = 0.0;

    switch (geometry_type)
      {
      case Cartesian :
	{
	for (unsigned int i = 0; i < dims; i++)
	  {
	  const double coord = logical_coords[i];
	  sq_dist += coord * coord;
	  }

	ret_val = 0;  // OK
	}
	break;

      case Polar :
	{
	const double dist = logical_coords[Polar_Radius_Index];

	sq_dist = dist * dist;

	ret_val = 0;  // OK
	}
	break;

      case Cylindrical :
	{
	const double xy_dist = logical_coords[Cylindrical_Radius_Index];

	sq_dist = xy_dist * xy_dist;

	if (dims > 1)
	  {
	  const double z = logical_coords[Cylindrical_zAxis_Index];

	  sq_dist += z * z;
	  }

	ret_val = 0;  // OK
	}
	break;

      case Spherical :
	{
	const double dist = logical_coords[Spherical_Radius_Index];

	sq_dist = dist * dist;
	
	ret_val = 0;  // OK
	}
	break;
      }

    origin_sqDistance = sq_dist;
    }

  return ret_val;
  }


int get_polar_axis_sqDist(const GeometryType geometry_type,
			  const std::vector<double> & logical_coords,
			  double & polar_axis_sqDistance)
  {
  int ret_val = -1;  // Error state

  const unsigned int dims = logical_coords.size();

  if ((dims >= 1) && (dims <= 3))
    {
    double sq_dist = 0.0;

    switch (geometry_type)
      {
      case Cartesian :
	{
	const double x = logical_coords[0];

	sq_dist = x * x;

	if (dims > 1)
	  {
	  const double y = logical_coords[1];
	  sq_dist += y * y;
	  }

	ret_val = 0;  // OK
	}
	break;

      case Polar :
	{
	const double dist = logical_coords[Polar_Radius_Index];

	sq_dist = dist * dist;

	ret_val = 0;  // OK
	}
	break;

      case Cylindrical :
	{
	const double xy_dist = logical_coords[Cylindrical_Radius_Index];

	sq_dist = xy_dist * xy_dist;

	ret_val = 0;  // OK
	}
	break;

      case Spherical :
	{
	const double dist = logical_coords[Spherical_Radius_Index];

	sq_dist = dist * dist;

	if (dims > 1)
	  {
	  const double polar = logical_coords[Spherical_Polar_Index];

	  const double sin_polar = sin(polar);

	  sq_dist *= sin_polar * sin_polar;
	  }
	
	ret_val = 0;  // OK
	}
	break;
      }

    polar_axis_sqDistance = sq_dist;
    }

  return ret_val;
  }




}  // End namespace Geometry
}  // End namespace QuickFlash
