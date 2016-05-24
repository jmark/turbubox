// C++ program file quickflash_block_cellvolume.cpp

/*
  By Nathan C. Hearn
     June 13, 2007

  Classes for computing cell volumes in various geometries.
*/


#include "quickflash_block_cellvolume.hpp"
#include <vector>
#include <math.h>
#include "quickflash_except.hpp"
#include "quickflash_geometry.hpp"


namespace QuickFlash
{
namespace Block
{  
namespace CellVolume
{

// VolumeBase object creator

const VolumeBase * new_volume_obj(const Geometry::GeometryType geom_type,
				  const std::vector<unsigned int> & block_dims,
				  const std::vector<double> & coord_deltas,
			   const std::vector<double> & min_cell_center_coords)
  {
  const VolumeBase * new_obj = 0;  // Default is null

  switch(geom_type)
    {
    case Geometry::Cartesian :
      new_obj = new CartesianVolume(coord_deltas);
      break;

    case Geometry::Polar :
      new_obj = new PolarVolume(block_dims, coord_deltas,
				min_cell_center_coords);
      break;

    case Geometry::Cylindrical :
      new_obj = new CylindricalVolume(block_dims, coord_deltas,
				      min_cell_center_coords);
      break;

    case Geometry::Spherical :
      new_obj = new SphericalVolume(block_dims, coord_deltas,
				    min_cell_center_coords);
      break;

    // DEFAULT CASE???
    }

  return new_obj;
  }


// Class CartesianVolume

void CartesianVolume::reset()
  { cell_vol = 0.0; }


void CartesianVolume::reset(const std::vector<double> & coord_deltas)
  {
  cell_vol = 1.0;

  const unsigned int dims = coord_deltas.size();

  for (unsigned int i = 0; i < dims; i++)
    cell_vol *= coord_deltas[i];
  }


void CartesianVolume::reset(const CartesianVolume & source)
  {
  if (&source != this)
    {
    cell_vol = source.cell_vol;
    }
  }


// Class PolarVolume

void PolarVolume::reset()
  {
  dims = 0;

  min_center_r = 0.0;

  delta_r = 0.0;
  delta_r_delta_azim = 0.0;

  azim_dims = 0;
  }


void PolarVolume::reset(const std::vector<unsigned int> & block_dims,
			const std::vector<double> & coord_deltas,
			const std::vector<double> & min_cell_center_coords)
  {
  reset();

  dims = block_dims.size();

  if ((dims > 2) || (dims < 1))
    throw Except("Incompatible dimensions for PolarVolume object", __FILE__,
		 __LINE__);

  if (dims == 2)
    {
    azim_dims = block_dims[Geometry::Polar_Azimuth_Index];

    if (azim_dims < 1)
      throw Except("Theta block dimensions must be at least one",
		   __FILE__, __LINE__);
    }
  else
    azim_dims = 1;

  if (coord_deltas.size() != dims)
    throw Except("Incompatible coordinate deltas size for PolarVolume object",
		 __FILE__, __LINE__);

  delta_r = coord_deltas[Geometry::Polar_Radius_Index];

  double delta_azim = 0.0;

  if (dims > 1)
    delta_azim = coord_deltas[Geometry::Polar_Azimuth_Index];
  else
    delta_azim = Geometry::two_pi;

  delta_r_delta_azim = delta_r * delta_azim;
  
  if (min_cell_center_coords.size() != dims)
    throw Except("Incompatible center coords size for PolarVolume object",
		 __FILE__, __LINE__);

  min_center_r = min_cell_center_coords[Geometry::Polar_Radius_Index];
  }


void PolarVolume::reset(const PolarVolume & source)
  {
  if (&source != this)
    {
    dims = source.dims;

    min_center_r = source.min_center_r;

    delta_r = source.delta_r;
    delta_r_delta_azim = source.delta_r_delta_azim;

    azim_dims = source.azim_dims;
    }
  }


double PolarVolume::cell_volume(const std::vector<unsigned int> & cell_index) 
  const
  {
  if (cell_index.size() != dims)
    throw Except("Incompatible cell_index vector for PolarVolume object",
		 __FILE__, __LINE__);

  const unsigned int r_index = cell_index[Geometry::Polar_Radius_Index];

  // ERROR CHECKING???

  return cell_volume_r(r_index);
  }


double PolarVolume:: cell_volume(const unsigned int cell_index) const
  {
  unsigned int r_index = cell_index;

  // ERROR CHECKING???

  if (dims == 2)
    r_index = cell_index / azim_dims;

  return cell_volume_r(r_index);
  }


// Class CylindricalVolume

void CylindricalVolume::reset()
  {
  dims = 0;

  min_center_r = 0.0;

  delta_r = 0.0;
  delta_r_delta_z_delta_azim = 0.0;

  z_dims_azim_dims = 0;
  }

 
void CylindricalVolume::reset(const std::vector<unsigned int> & block_dims,
			      const std::vector<double> & coord_deltas,
			  const std::vector<double> & min_cell_center_coords)
  {
  reset();

  dims = block_dims.size();

  if ((dims > 3) || (dims < 1))
    throw Except("Incompatible dimensions for CylindricalVolume object",
		 __FILE__, __LINE__);

  z_dims_azim_dims = 1;

  if (dims > 1)
    {
    const unsigned int z_dims = block_dims[Geometry::Cylindrical_zAxis_Index];

    if (z_dims < 1)
      throw Except("z-axis block dimensions must be at least one",
		   __FILE__, __LINE__);

    z_dims_azim_dims *= z_dims;
    }

  if (dims > 2)
    {
    const unsigned int azim_dims 
      = block_dims[Geometry::Cylindrical_Azimuth_Index];

    if (azim_dims < 1)
      throw Except("Azimuthal block dimensions must be at least one");

    z_dims_azim_dims *= azim_dims;
    }

  if (coord_deltas.size() != dims)
    throw Except("Incompatible coord deltas size for CylindricalVolume object",
		 __FILE__, __LINE__);

  delta_r = coord_deltas[Geometry::Cylindrical_Radius_Index];

  double delta_z = 1.0;

  if (dims > 1)
    delta_z = coord_deltas[Geometry::Cylindrical_zAxis_Index];

  double delta_azim = 0.0;

  if (dims > 2)
    delta_azim = coord_deltas[Geometry::Cylindrical_Azimuth_Index];
  else
    delta_azim = Geometry::two_pi;

  delta_r_delta_z_delta_azim = delta_r * delta_z * delta_azim;
  
  if (min_cell_center_coords.size() != dims)
    throw Except("Incompatible coords size for CylindricalVolume object",
		 __FILE__, __LINE__);

  min_center_r = min_cell_center_coords[Geometry::Cylindrical_Radius_Index];
  }


void CylindricalVolume::reset(const CylindricalVolume & source)
  {
  if (&source != this)
    {
    dims = source.dims;

    min_center_r = source.min_center_r;

    delta_r = source.delta_r;
    delta_r_delta_z_delta_azim = source.delta_r_delta_z_delta_azim;

    z_dims_azim_dims = source.z_dims_azim_dims;
    }
  }


double CylindricalVolume::cell_volume(const std::vector<unsigned int> & 
				      cell_index) const
  {
  if (cell_index.size() != dims)
    throw Except("Incompatible cell_index vector for CylindricalVolume object",
		 __FILE__, __LINE__);

  const unsigned int r_index = cell_index[Geometry::Cylindrical_Radius_Index];

  // ERROR CHECKING???

  return cell_volume_r(r_index);
  }


double CylindricalVolume:: cell_volume(const unsigned int cell_index) const
  {
  unsigned int r_index = cell_index;

  // ERROR CHECKING???

  if (dims > 1)
    r_index = cell_index / z_dims_azim_dims;

  return cell_volume_r(r_index);
  }


// Class SphericalVolume

const double SphericalVolume::one_twelfth = 1.0 / 12.0;


void SphericalVolume::reset()
  {
  dims = 0;

  min_center_r = 0.0;

  delta_r = 0.0;
  delta_r_squared_over_twelve = 0.0;

  delta_azim = 0.0;

  polar_cell_volcomp.clear();

  azim_dims = 0;
  polar_dims = 0;
  }


void SphericalVolume::reset(const std::vector<unsigned int> & block_dims,
			    const std::vector<double> & coord_deltas,
			    const std::vector<double> & min_cell_center_coords)
  {
  reset();

  dims = block_dims.size();

  if ((dims > 3) || (dims < 1))
    throw Except("Incompatible dimensions for SphericalVolume object", 
		 __FILE__, __LINE__);

  polar_dims = 1;
  azim_dims = 1;

  if (dims > 1)
    {
    polar_dims = block_dims[Geometry::Spherical_Polar_Index];

    if (polar_dims < 1)
      throw Except("Polar block dimensions must be at least one",
		   __FILE__, __LINE__);
    }

  if (dims > 2)
    {
    azim_dims = block_dims[Geometry::Spherical_Azimuth_Index];

    if (azim_dims < 1)
      throw Except("Azimuthal block dimensions must be at least one");
    }

  if (min_cell_center_coords.size() != dims)
    throw Except("Incompatible coords size for CylindricalVolume object",
		 __FILE__, __LINE__);

  min_center_r = min_cell_center_coords[Geometry::Spherical_Radius_Index];

  if (coord_deltas.size() != dims)
    throw Except("Incompatible coord deltas size for SphericalVolume object",
		 __FILE__, __LINE__);

  delta_r = coord_deltas[Geometry::Spherical_Radius_Index];

  delta_r_squared_over_twelve = one_twelfth * (delta_r * delta_r);

  polar_cell_volcomp.resize(polar_dims);

  if (dims > 1)
    {
    const double min_center_polar 
      = min_cell_center_coords[Geometry::Spherical_Polar_Index];

    const double delta_polar = coord_deltas[Geometry::Spherical_Polar_Index];

    double current_cell_min_polar = min_center_polar - (delta_polar * 0.5);

    // THE NUMBER OF COS CALLS COULD BE REDUCED BY KEEPING THE PREVIOUS VALUE

    for (unsigned int idx = 0; idx < polar_dims; idx++)
      {
      const double current_cell_max_polar 
	= current_cell_min_polar + delta_polar;

      const double polar_vol_component 
	= cos(current_cell_min_polar) - cos(current_cell_max_polar);

      polar_cell_volcomp[idx] = polar_vol_component;

      current_cell_min_polar = current_cell_max_polar;
      }
    }
  else
    polar_cell_volcomp[0] = 2.0;  // Full sphere

  if (dims > 2)
    delta_azim = coord_deltas[Geometry::Spherical_Azimuth_Index];
  else
    delta_azim = Geometry::two_pi;
  }


void SphericalVolume::reset(const SphericalVolume & source)
  {
  if (&source != this)
    {
    dims = source.dims;

    min_center_r = source.min_center_r;

    delta_r = source.delta_r;
    delta_r_squared_over_twelve = source.delta_r_squared_over_twelve;

    delta_azim = source.delta_azim;

    polar_cell_volcomp = source.polar_cell_volcomp;

    polar_dims = source.polar_dims;
    azim_dims = source.azim_dims;
    }
  }


double SphericalVolume::cell_volume(const std::vector<unsigned int> & 
				    cell_index) const
  {
  if (cell_index.size() != dims)
    throw Except("Incompatible cell_index vector for SphericalVolume object",
		 __FILE__, __LINE__);

  const unsigned int r_index = cell_index[Geometry::Spherical_Radius_Index];

  unsigned int polar_index = 0;

  if (dims > 1)
    polar_index = cell_index[Geometry::Spherical_Polar_Index];

  // ERROR CHECKING???

  return cell_volume_r_polar(r_index, polar_index);
  }


double SphericalVolume::cell_volume(const unsigned int cell_index) const
  {
  const unsigned int no_azim_index = cell_index / azim_dims;

  const unsigned int polar_index = no_azim_index % polar_dims;

  const unsigned int r_index = no_azim_index / polar_dims;

  // ERROR CHECKING???

  return cell_volume_r_polar(r_index, polar_index);
  }


double SphericalVolume::cell_volume_r_polar(const unsigned int r_index,
					    const unsigned int polar_index) 
  const
  {
  const double r = min_center_r + (static_cast<double>(r_index) * delta_r);

  const double polar_vol_comp = polar_cell_volcomp[polar_index];

  const double r_vol_comp = 
    ((r * r) + delta_r_squared_over_twelve) * delta_r;

  return r_vol_comp * polar_vol_comp * delta_azim;
  }


}  // End namespace CellVolume
}  // End namespace Block
}  // End namespace QuickFlash
