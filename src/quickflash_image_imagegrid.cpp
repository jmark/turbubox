// C++ program file quickflash_image_imagegrid.cpp

/*
  By Nathan C. Hearn
     February 11, 2008

  General floating point image class.
*/


#include "quickflash_image_imagegrid.hpp"
#include <vector>
#include "quickflash_except.hpp"
#include "quickflash_block_blockdata.hpp"
#include "quickflash_color_color.hpp"
#include "quickflash_color_colormap.hpp"


namespace QuickFlash
{
namespace Image
{

// Class ImageGrid

ImageGrid::ImageGrid() : image() { }


ImageGrid::ImageGrid(const std::vector<unsigned int> & pixel_count) :
  image()
  { reset(pixel_count); }

ImageGrid::ImageGrid(const std::vector<unsigned int> & pixel_count, 
		     const Color::Color & color) :
  image()
  { reset(pixel_count, color); }

ImageGrid::ImageGrid(const Block::BlockData<double> & data, 
		     const Color::Colormap & colormap) :
  image()
  { reset(data, colormap); }

ImageGrid::ImageGrid(const ImageGrid & source) :
  image()
  { reset(source); }


void ImageGrid::reset()
  { image.reset(); }


void ImageGrid::reset(const std::vector<unsigned int> & pixel_count)
  { 
  const Color::Color default_color;

  reset(pixel_count, default_color);
  }


void ImageGrid::reset(const std::vector<unsigned int> & pixel_count, 
		      const Color::Color & color)
  {
  image.reset(pixel_count);

  reset_pixels(color);
  }


void ImageGrid::reset(const Block::BlockData<double> & data, 
		      const Color::Colormap & colormap)
  {
  const unsigned int num_pixels = data.get_num_cells();

  if (num_pixels < 1)
    reset();
  else
    {
    const unsigned int num_dims = data.get_space_dims();

    if (num_dims != 2)
      throw Except("Incompatible dimensions for ImageGrid", __FILE__,
		   __LINE__);

    const std::vector<unsigned int> & image_dims = data.get_block_dims();

    image.reset(image_dims);

    // Fill in the color values

    for (unsigned int index = 0; index < num_pixels; index++)
      {
      const double cell_data = data[index];

      colormap.get_color(cell_data, image[index]);
      }
    }
  }


void ImageGrid::reset(const ImageGrid & source)
  {
  if (&source != this)
    {
    image = source.image;
    }
  }


void ImageGrid::reset_pixels(const Color::Color & color)
  {
  const unsigned int num_pixels = image.get_num_cells();

  for (unsigned int index = 0; index < num_pixels; index++)
    image[index] = color;
  }


void ImageGrid::set_alpha(const double alpha_value)
  {
  const unsigned int num_pixels = image.get_num_cells();

  for (unsigned int index = 0; index < num_pixels; index++)
    image[index].set_alpha(alpha_value);
  }


void ImageGrid::add_image(const ImageGrid & other_image, 
			  const double current_scale, const double other_scale)
  {
  const std::vector<int> image_offset(2, 0);

  add_image(other_image, image_offset, current_scale, other_scale);
  }


void ImageGrid::add_image(const ImageGrid & other_image, 
			  const std::vector<int> & other_image_offset,
			  const double current_scale, const double other_scale)
  {
  const std::vector<unsigned int> & other_image_dims 
    = other_image.get_image_dims();

  std::vector<unsigned int> offset_onto_current_image;
  std::vector<unsigned int> other_image_start;
  std::vector<unsigned int> other_image_pixel_count;

  get_overlay_bounds(other_image_dims, other_image_offset,
		     offset_onto_current_image, other_image_start,
		     other_image_pixel_count);

  const unsigned int num_i = other_image_pixel_count[0];
  const unsigned int num_j = other_image_pixel_count[1];

  std::vector<unsigned int> current_image_coords(2);
  std::vector<unsigned int> other_image_coords(2);

  for (unsigned int i = 0; i < num_i; i++)
    {
    current_image_coords[0] = i + offset_onto_current_image[0];
    other_image_coords[0] = i + other_image_start[0];

    for (unsigned int j = 0; j < num_j; j++)
      {
      current_image_coords[1] = j + offset_onto_current_image[1];
      other_image_coords[1] = j + other_image_start[1];

      const Color::Color & other_color = other_image[other_image_coords];

      image[current_image_coords].blend_color(other_color, current_scale, 
					      other_scale);
      }
    }
  }


void ImageGrid::add_image(const Block::BlockData<double> & other_data, 
			  const Color::Colormap & colormap,
			  const double current_scale, const double other_scale)
  {
  const std::vector<int> image_offset(2, 0);

  add_image(other_data, colormap, image_offset, current_scale, other_scale);
  }


void ImageGrid::add_image(const Block::BlockData<double> 
			  & other_data, 
			  const Color::Colormap & colormap,
			  const std::vector<int> & other_image_offset,
			  const double current_scale, const double other_scale)
  {
  const unsigned int num_dims = other_data.get_space_dims();

  if (num_dims != 2)
    throw Except("Incompatible dimensions for ImageGrid", __FILE__, __LINE__);

  const std::vector<unsigned int> & other_image_dims 
    = other_data.get_block_dims();

  std::vector<unsigned int> offset_onto_current_image;
  std::vector<unsigned int> other_image_start;
  std::vector<unsigned int> other_image_pixel_count;

  get_overlay_bounds(other_image_dims, other_image_offset,
		     offset_onto_current_image, other_image_start,
		     other_image_pixel_count);

  const unsigned int num_i = other_image_pixel_count[0];
  const unsigned int num_j = other_image_pixel_count[1];

  std::vector<unsigned int> current_image_coords(2);
  std::vector<unsigned int> other_image_coords(2);

  Color::Color other_color;

  for (unsigned int i = 0; i < num_i; i++)
    {
    current_image_coords[0] = i + offset_onto_current_image[0];
    other_image_coords[0] = i + other_image_start[0];

    for (unsigned int j = 0; j < num_j; j++)
      {
      current_image_coords[1] = j + offset_onto_current_image[1];
      other_image_coords[1] = j + other_image_start[1];

      const double other_value = other_data[other_image_coords];

      colormap.get_color(other_value, other_color);

      image[current_image_coords].blend_color(other_color, current_scale, 
					      other_scale);
      }
    }
  }


void ImageGrid::overlay_image(const ImageGrid & other_image, 
			      const bool add_below)
  {
  const std::vector<int> image_offset(2, 0);

  overlay_image(other_image, image_offset, add_below);
  }


void ImageGrid::overlay_image(const ImageGrid & other_image, 
			      const std::vector<int> & other_image_offset,
			      const bool add_below)
  {
  const std::vector<unsigned int> & other_image_dims 
    = other_image.get_image_dims();

  std::vector<unsigned int> offset_onto_current_image;
  std::vector<unsigned int> other_image_start;
  std::vector<unsigned int> other_image_pixel_count;

  get_overlay_bounds(other_image_dims, other_image_offset,
		     offset_onto_current_image, other_image_start,
		     other_image_pixel_count);

  const unsigned int num_i = other_image_pixel_count[0];
  const unsigned int num_j = other_image_pixel_count[1];

  std::vector<unsigned int> current_image_coords(2);
  std::vector<unsigned int> other_image_coords(2);

  for (unsigned int i = 0; i < num_i; i++)
    {
    current_image_coords[0] = i + offset_onto_current_image[0];
    other_image_coords[0] = i + other_image_start[0];

    for (unsigned int j = 0; j < num_j; j++)
      {
      current_image_coords[1] = j + offset_onto_current_image[1];
      other_image_coords[1] = j + other_image_start[1];

      const Color::Color & other_color = other_image[other_image_coords];

      image[current_image_coords].overlay_color(other_color, add_below);
      }
    }
  }


void ImageGrid::overlay_image(const Block::BlockData<double> 
			      & other_data, 
			      const Color::Colormap & colormap,
			      const bool add_below)
  {
  const std::vector<int> image_offset(2, 0);

  overlay_image(other_data, colormap, image_offset, add_below);
  }


void ImageGrid::overlay_image(const Block::BlockData<double> 
			      & other_data, 
			      const Color::Colormap & colormap,
			      const std::vector<int> & other_image_offset,
			      const bool add_below)
  {
  const unsigned int num_dims = other_data.get_space_dims();

  if (num_dims != 2)
    throw Except("Incompatible dimensions for ImageGrid", __FILE__, __LINE__);

  const std::vector<unsigned int> & other_image_dims 
    = other_data.get_block_dims();

  std::vector<unsigned int> offset_onto_current_image;
  std::vector<unsigned int> other_image_start;
  std::vector<unsigned int> other_image_pixel_count;

  get_overlay_bounds(other_image_dims, other_image_offset,
		     offset_onto_current_image, other_image_start,
		     other_image_pixel_count);

  const unsigned int num_i = other_image_pixel_count[0];
  const unsigned int num_j = other_image_pixel_count[1];

  std::vector<unsigned int> current_image_coords(2);
  std::vector<unsigned int> other_image_coords(2);

  Color::Color other_color;

  for (unsigned int i = 0; i < num_i; i++)
    {
    current_image_coords[0] = i + offset_onto_current_image[0];
    other_image_coords[0] = i + other_image_start[0];

    for (unsigned int j = 0; j < num_j; j++)
      {
      current_image_coords[1] = j + offset_onto_current_image[1];
      other_image_coords[1] = j + other_image_start[1];

      const double other_value = other_data[other_image_coords];

      colormap.get_color(other_value, other_color);

      image[current_image_coords].overlay_color(other_color, add_below);
      }
    }
  }


void ImageGrid::get_overlay_bounds(const std::vector<unsigned int> 
				     & other_image_dims,
				   const std::vector<int> & other_image_offset,
				   std::vector<unsigned int> 
				     & offset_onto_current_image,
				   std::vector<unsigned int> 
				     & other_image_start,
				   std::vector<unsigned int> 
				     & other_image_pixel_count) const
  {
  if (other_image_offset.size() != 2)
    throw Except("Incompatible offset vector", __FILE__, __LINE__);

  const std::vector<unsigned int> & current_image_dims 
    = image.get_block_dims();

  if (current_image_dims.size() != 2)
    throw Except("Current image not initialized", __FILE__, __LINE__);

  if (other_image_dims.size() != 2)
    throw Except("Overlay image not initialized", __FILE__, __LINE__);

  // Determine the area of overlap

  offset_onto_current_image.resize(2);
  other_image_start.resize(2);
  other_image_pixel_count.resize(2);

  for (unsigned int i = 0; i < 2; i++)
    {
    const unsigned int current_image_axis_size = current_image_dims[i];

    const unsigned int other_image_axis_size = other_image_dims[i];

    const int other_image_axis_offset = other_image_offset[i];

    unsigned int current_image_axis_offset = 0;

    unsigned int other_image_axis_start = 0;

    if (other_image_axis_offset < 0)
      other_image_axis_start 
	= static_cast<unsigned int>(-other_image_axis_offset);
    else
      current_image_axis_offset 
	= static_cast<unsigned int>(other_image_axis_offset);

    unsigned int other_image_axis_count 
      = other_image_axis_size - other_image_axis_start;

    if (current_image_axis_size 
	< (other_image_axis_count + current_image_axis_offset))
      {
      other_image_axis_count 
	= current_image_axis_size - current_image_axis_offset;
      }

    offset_onto_current_image[i] = current_image_axis_offset;
    other_image_start[i] = other_image_axis_start;
    other_image_pixel_count[i] = other_image_axis_count;
    }
  }


}  // End namespace Image
}  // End namespace QuickFlash
