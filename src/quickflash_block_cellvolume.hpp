// C++ header file quickflash_block_cellvolume.hpp

/*
  By Nathan C. Hearn
     June 11, 2007

  Classes for computing cell volumes in various geometries.
*/


#ifndef QUICKFLASH_BLOCK_CELLVOLUME_HPP
#define QUICKFLASH_BLOCK_CELLVOLUME_HPP


#include <vector>
#include "quickflash_geometry.hpp"


namespace QuickFlash
{
namespace Block
{  
namespace CellVolume
{

// Cell volume base class prototype

class VolumeBase;


// VolumeBase object creator

const VolumeBase * new_volume_obj(const Geometry::GeometryType geom_type,
				  const std::vector<unsigned int> & block_dims,
				  const std::vector<double> & coord_deltas,
			   const std::vector<double> & min_cell_center_coords);


// Cell volume base class

class VolumeBase
  {
  public:
    VolumeBase() { }
    virtual ~VolumeBase() { }

    virtual VolumeBase * clone() const = 0;

    virtual double cell_volume(const std::vector<unsigned int> & cell_index) 
      const = 0;

    virtual double cell_volume(const unsigned int cell_index) const = 0;
  };



// Class CartesianVolume

class CartesianVolume : public VolumeBase
  {
  public:
    CartesianVolume() : VolumeBase(), cell_vol(0.0) { }

    CartesianVolume(const std::vector<double> & coord_deltas) :
      VolumeBase(), cell_vol(0.0)
      { reset(coord_deltas); }

    CartesianVolume(const CartesianVolume & source) :
      VolumeBase(), cell_vol(0.0)
      { reset(source); }

    virtual ~CartesianVolume() { }

    CartesianVolume & operator=(const CartesianVolume & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const std::vector<double> & coord_deltas);

    void reset(const CartesianVolume & source);

    virtual VolumeBase * clone() const
      { return new CartesianVolume(*this); }

    virtual double cell_volume(const std::vector<unsigned int> &) const
      { return cell_vol; }

    virtual double cell_volume(const unsigned int) const
      { return cell_vol; }

  private:
    double cell_vol;
  };


// Class PolarVolume

class PolarVolume : public VolumeBase
  {
  public:
    PolarVolume() :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_delta_azim(0.0), azim_dims(0)
      { }

    PolarVolume(const std::vector<unsigned int> & block_dims,
		const std::vector<double> & coord_deltas,
		const std::vector<double> & min_cell_center_coords) :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_delta_azim(0.0), azim_dims(0)
      { reset(block_dims, coord_deltas, min_cell_center_coords); }

    PolarVolume(const PolarVolume & source) :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_delta_azim(0.0), azim_dims(0)
      { reset(source); }

    virtual ~PolarVolume() { }

    PolarVolume & operator=(const PolarVolume & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const std::vector<unsigned int> & block_dims,
	       const std::vector<double> & coord_deltas,
	       const std::vector<double> & min_cell_center_coords);

    void reset(const PolarVolume & source);

    virtual VolumeBase * clone() const
      { return new PolarVolume(*this); }

    virtual double cell_volume(const std::vector<unsigned int> & cell_index) 
      const;

    virtual double cell_volume(const unsigned int cell_index) const;

  private:
    double cell_volume_r(const unsigned int r_index) const
      {  
      const double r = min_center_r + (static_cast<double>(r_index) * delta_r);
      return r * delta_r_delta_azim;
      }

  private:
    unsigned int dims;

    double min_center_r;

    double delta_r;
    double delta_r_delta_azim;

    unsigned int azim_dims;
  };


class CylindricalVolume : public VolumeBase
  {
  public:
    CylindricalVolume() :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_delta_z_delta_azim(0.0), z_dims_azim_dims(0)
      { }

    CylindricalVolume(const std::vector<unsigned int> & block_dims,
		const std::vector<double> & coord_deltas,
		const std::vector<double> & min_cell_center_coords) :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_delta_z_delta_azim(0.0), z_dims_azim_dims(0)
      { reset(block_dims, coord_deltas, min_cell_center_coords); }

    CylindricalVolume(const CylindricalVolume & source) :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_delta_z_delta_azim(0.0), z_dims_azim_dims(0)
      { reset(source); }

    virtual ~CylindricalVolume() { }

    CylindricalVolume & operator=(const CylindricalVolume & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const std::vector<unsigned int> & block_dims,
	       const std::vector<double> & coord_deltas,
	       const std::vector<double> & min_cell_center_coords);

    void reset(const CylindricalVolume & source);

    virtual VolumeBase * clone() const
      { return new CylindricalVolume(*this); }

    virtual double cell_volume(const std::vector<unsigned int> & cell_index) 
      const;

    virtual double cell_volume(const unsigned int cell_index) const;

  private:
    double cell_volume_r(const unsigned int r_index) const
      {  
      const double r = min_center_r + (static_cast<double>(r_index) * delta_r);
      return r * delta_r_delta_z_delta_azim;
      }
				       
  private:
    unsigned int dims;

    double min_center_r;

    double delta_r;
    double delta_r_delta_z_delta_azim;

    unsigned int z_dims_azim_dims;
  };


class SphericalVolume : public VolumeBase
  {
  public:
    SphericalVolume() :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_squared_over_twelve(0.0), delta_azim(0.0), polar_cell_volcomp(),
      polar_dims(0), azim_dims(0)
      { }

    SphericalVolume(const std::vector<unsigned int> & block_dims,
		    const std::vector<double> & coord_deltas,
		    const std::vector<double> & min_cell_center_coords) :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_squared_over_twelve(0.0), delta_azim(0.0), polar_cell_volcomp(),
      polar_dims(0), azim_dims(0)
      { reset(block_dims, coord_deltas, min_cell_center_coords); }

    SphericalVolume(const SphericalVolume & source) :
      VolumeBase(), dims(0), min_center_r(0.0), delta_r(0.0), 
      delta_r_squared_over_twelve(0.0), delta_azim(0.0), polar_cell_volcomp(),
      polar_dims(0), azim_dims(0)
      { reset(source); }

    virtual ~SphericalVolume() { }

    SphericalVolume & operator=(const SphericalVolume & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const std::vector<unsigned int> & block_dims,
	       const std::vector<double> & coord_deltas,
	       const std::vector<double> & min_cell_center_coords);

    void reset(const SphericalVolume & source);

    virtual VolumeBase * clone() const
      { return new SphericalVolume(*this); }

    virtual double cell_volume(const std::vector<unsigned int> & cell_index) 
      const;

    virtual double cell_volume(const unsigned int cell_index) const;

  private:
    double cell_volume_r_polar(const unsigned int r_index,
			       const unsigned int polar_index=0) const;

  private:
    static const double one_twelfth;

    unsigned int dims;

    double min_center_r;

    double delta_r;
    double delta_r_squared_over_twelve;

    double delta_azim;

    std::vector<double> polar_cell_volcomp;

    unsigned int polar_dims;
    unsigned int azim_dims;
  };


}  // End namespace CellVolume
}  // End namespace Block
}  // End namespace QuickFlash


#endif  // QUICKFLASH_BLOCK_CELLVOLUME_HPP
