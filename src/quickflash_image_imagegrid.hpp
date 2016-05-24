// C++ header file quickflash_image_imagegrid.hpp

/*
  By Nathan C. Hearn
     February 11, 2008

  General floating point image class.
*/


#ifndef QUICKFLASH_IMAGE_IMAGEGRID_HPP
#define QUICKFLASH_IMAGE_IMAGEGRID_HPP


#include <vector>
#include "quickflash_block_blockdata.hpp"
#include "quickflash_color_color.hpp"
#include "quickflash_color_colormap.hpp"


namespace QuickFlash
{
namespace Image
{

class ImageGrid
  {
  public :

    ImageGrid();

    ImageGrid(const std::vector<unsigned int> & pixel_count);

    ImageGrid(const std::vector<unsigned int> & pixel_count, 
	      const Color::Color & color);

    ImageGrid(const Block::BlockData<double> & data, 
	      const Color::Colormap & colormap);

    ImageGrid(const ImageGrid & source);

    ImageGrid & operator=(const ImageGrid & source)
      {
      reset(source);
      return *this;
      }

    Color::Color & operator[](const unsigned int pixel_index)
      { return image[pixel_index]; }

    const Color::Color & operator[](const unsigned int pixel_index) const
      { return image[pixel_index]; }

    Color::Color & operator[](const std::vector<unsigned int> & pixel_index)
      { return image[pixel_index]; }

    const Color::Color & operator[](const std::vector<unsigned int> 
				    & pixel_index) const
      { return image[pixel_index]; }

    void reset();

    void reset(const std::vector<unsigned int> & pixel_count);
    void reset(const std::vector<unsigned int> & pixel_count, 
	       const Color::Color & color);

    void reset(const Block::BlockData<double> & data, 
	       const Color::Colormap & colormap);

    void reset(const ImageGrid & source);

    void reset_pixels(const Color::Color & color);

    unsigned int get_num_pixels() const { return image.get_num_cells(); }

    const std::vector<unsigned int> & get_image_dims() const 
      { return image.get_block_dims(); }

    void set_alpha(const double alpha_value);

    void add_image(const ImageGrid & other_image, 
		   const double current_scale=1.0,
		   const double other_scale=1.0);

    void add_image(const ImageGrid & other_image, 
		   const std::vector<int> & other_image_offset,
		   const double current_scale=1.0, 
		   const double other_scale=1.0);

    void add_image(const Block::BlockData<double> & other_data, 
		   const Color::Colormap & colormap,
		   const double current_scale=1.0, 
		   const double other_scale=1.0);

    void add_image(const Block::BlockData<double> & other_data, 
		   const Color::Colormap & colormap,
		   const std::vector<int> & other_image_offset,
		   const double current_scale=1.0,
		   const double other_scale=1.0);

    void overlay_image(const ImageGrid & other_image, 
		       const bool add_below=false);

    void overlay_image(const ImageGrid & other_image, 
		       const std::vector<int> & other_image_offset,
		       const bool add_below=false);

    void overlay_image(const Block::BlockData<double> & other_data, 
		       const Color::Colormap & colormap,
		       const bool add_below=false);

    void overlay_image(const Block::BlockData<double> & other_data, 
		       const Color::Colormap & colormap,
		       const std::vector<int> & other_image_offset,
		       const bool add_below=false);

  private :

    void get_overlay_bounds(const std::vector<unsigned int> & other_image_dims,
			    const std::vector<int> & other_image_offset,
			    std::vector<unsigned int> 
			      & offset_onto_current_image,
			    std::vector<unsigned int> & other_image_start,
			    std::vector<unsigned int> 
			      & other_image_pixel_count) const;

  protected :

    Block::BlockData<Color::Color> image;
  };


}  // End namespace Image
}  // End namespace QuickFlash


#endif  // QUICKFLASH_IMAGE_IMAGEGRID_HPP
