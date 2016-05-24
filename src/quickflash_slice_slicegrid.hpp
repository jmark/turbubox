// C++ header file quickflash_slice_slicegrid.hpp

/*
  By Nathan C. Hearn
     February 11, 2008

  Grid for composing 2D slices.
*/


#ifndef QUICKFLASH_SLICE_SLICEGRID_HPP
#define QUICKFLASH_SLICE_SLICEGRID_HPP

#include <vector>
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_block_blockdata.hpp"


namespace QuickFlash
{
namespace Slice
{

class SliceGrid
  {
  public :
    SliceGrid() : 
      space_dims(0), grid_info(), grid_data(), grid_center(), 
      slice_width(), orthog_min(0.0), orthog_width(0.0), 
      roll(0.0), polar(0.0), azim(0.0), trans_mat() 
    { }

    SliceGrid(const std::vector<double> & center_coords, 
	      const std::vector<double> & slice_grid_width, 
	      const std::vector<unsigned int> & pixel_count, 
	      const double azimuthal_angle=0.0, const double polar_angle=0.0, 
	      const double roll_angle=0.0, const double orthogonal_width=0.0) :
      space_dims(0), grid_info(), grid_data(), grid_center(), 
      slice_width(), orthog_min(0.0), orthog_width(0.0), 
      roll(0.0), polar(0.0), azim(0.0), trans_mat()
      { 
      reset(center_coords, slice_grid_width, pixel_count, azimuthal_angle, 
	    polar_angle, roll_angle, orthogonal_width);
      }

    SliceGrid(const SliceGrid & source) :
      space_dims(0), grid_info(), grid_data(), grid_center(), 
      slice_width(), orthog_min(0.0), orthog_width(0.0), 
      roll(0.0), polar(0.0), azim(0.0), trans_mat()
      { reset(source); }

    ~SliceGrid() { }

    SliceGrid & operator=(const SliceGrid & source)
      {
      reset(source);
      return *this;
      }

    double & operator[](const unsigned int pixel_index)
      { return grid_data[pixel_index]; }

    double operator[](const unsigned int pixel_index) const
      { return grid_data[pixel_index]; }

    double & operator[](const std::vector<unsigned int> & pixel_index)
      { return grid_data[pixel_index]; }

    double operator[](const std::vector<unsigned int> & pixel_index) const
      { return grid_data[pixel_index]; }

    void reset();

    void reset(const std::vector<double> & center_coords, 
	       const std::vector<double> & slice_grid_width, 
	       const std::vector<unsigned int> & pixel_count, 
	       const double azimuthal_angle=0.0, const double polar_angle=0.0, 
	       const double roll_angle=0.0, const double orthogonal_width=0.0);

    void reset(const SliceGrid & source);

    void reset_metadata(const SliceGrid & source);

    unsigned int get_space_dims() const { return space_dims; }

    unsigned int get_num_pixels() const { return grid_info.get_num_cells(); }

    const std::vector<double> & get_grid_center() const { return grid_center; }

    const std::vector<unsigned int> & get_grid_dims() const 
      { return grid_data.get_block_dims(); }

    const std::vector<double> & get_slice_grid_width() const
      { return slice_width; }

    double get_orthogonal_width() const { return orthog_width; }

    double get_roll_angle() const { return roll; }
    double get_polar_angle() const { return polar; }
    double get_azimuthal_angle() const { return azim; }
    
    void get_minmax_values(double & min_value, double & max_value) const
      { grid_data.get_minmax_values(min_value, max_value); }

    void get_pixel_center(const unsigned int pixel_index,
			  std::vector<double> & space_coords) const;

    void get_pixel_center(const std::vector<unsigned int> & pixel_index,
			  std::vector<double> & space_coords) const;

    void get_pixel_subsamples(const unsigned int pixel_index, 
			      const std::vector<unsigned int> & sample_dims,
			      std::vector< std::vector<double> > 
			      & sample_space_coords) const;

    void get_pixel_subsamples(const std::vector<unsigned int> & pixel_index, 
			      const std::vector<unsigned int> & sample_dims,
			      std::vector< std::vector<double> > 
			      & sample_space_coords) const;

    void get_symmetry_pair(const unsigned int symmetry_axis, 
			   const unsigned int pixel_index,
			   double & symm_sum, double & symm_diff) const;

    void get_symmetry_pair(const unsigned int symmetry_axis, 
			   const std::vector<unsigned int> & pixel_index,
			   double & symm_sum, double & symm_diff) const;

    void set_symmetry_sum(const unsigned int symmetry_axis, 
			  const SliceGrid & source);

    void set_symmetry_diff(const unsigned int symmetry_axis,
			   const SliceGrid & source);

    const Block::BlockInfo & get_grid_info() const 
      { return grid_info; }

    const Block::BlockData<double> & get_grid_data() const 
      { return grid_data; }

  private :
    void set_transform_matrix();
    void transform_pos(const std::vector<double> & grid_coords,
		       std::vector<double> & space_coords) const;

  private :
    unsigned int space_dims;

    Block::BlockInfo grid_info;
    Block::BlockData<double> grid_data;

    std::vector<double> grid_center;

    std::vector<double> slice_width;

    double orthog_min;
    double orthog_width;

    double roll;
    double polar;
    double azim;

    std::vector< std::vector<double> > trans_mat;
  };


}  // End namespace Slice
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SLICE_SLICEGRID_HPP
