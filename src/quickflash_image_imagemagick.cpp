// C++ program file quickflash_image_imagemagick.cpp

/*
  By Nathan C. Hearn
     February 12, 2008

  Colormap legend routines by John Norris.

  Routines for converting ImageGrid object data into bitmapped images.

  Requires ImageMagick's Magick++ bindings.
*/


#ifdef USE_MAGICK


#include "quickflash_image_imagemagick.hpp"
#include <Magick++.h>
#include <string>
#include <algorithm>
#include <fstream>
#include <cmath>
#include "quickflash_except.hpp"
#include "quickflash_image_imagegrid.hpp"
#include "quickflash_color_color.hpp"
#include "quickflash_color_colormap.hpp"
#include "quickflash_utils_file.hpp"


namespace QuickFlash
{
namespace Image
{
namespace ImageMagick
{

int write_image(const std::string & filename, const ImageGrid & data,
		const bool transpose_image)
  {
  int ret_val = -1;  // Error state

  if (filename.size() < 1)
    throw Except("Empty filename", __FILE__, __LINE__);

  const unsigned int num_pixels = data.get_num_pixels();

  if (num_pixels < 1)
    throw Except("No data for image", __FILE__, __LINE__);

  const std::vector<unsigned int> & image_dims = data.get_image_dims();

  // Determine which axis is which based on transpose_image flag

  const unsigned int horizontal_pixels 
    = transpose_image ? image_dims[1] : image_dims[0];

  const unsigned int vertical_pixels 
    = transpose_image ? image_dims[0] : image_dims[1];

  const Magick::Geometry geometry(horizontal_pixels, vertical_pixels);

  Magick::ColorRGB magick_color(0.0, 0.0, 0.0);

  Magick::Image image(geometry, magick_color);

  image.type(Magick::TrueColorMatteType);
  image.depth(8);

  Magick::Pixels pixel_cache(image);

  Magick::PixelPacket * pixels = pixel_cache.get(0, 0, horizontal_pixels, 
						 vertical_pixels);

  // NOTE: Image is stored beginning with upper-left corner, proceding in rows

  std::vector<unsigned int> index(2);

  for (unsigned int j = vertical_pixels; j > 0; j--)
    {
    if (transpose_image)
      index[0] = j - 1;
    else
      index[1] = j - 1;

    for (unsigned int i = 0; i < horizontal_pixels; i++)
      {
      if (transpose_image)
	index[1] = i;
      else
	index[0] = i;

      const Color::Color & color = data[index];

      magick_color.red(color.get_red());
      magick_color.green(color.get_green());
      magick_color.blue(color.get_blue());

      // ImageMagick (prior to 6.3.0) defines alpha as level of transparency
      // ... Magick++ seems to follow this definition for now ...

      const double transparency = 1.0 - color.get_alpha();

      magick_color.alpha(transparency);

      *(pixels) = magick_color;

      pixels++;
      }
    }

  image.syncPixels();

  image.write(filename);

  ret_val = 0;  // OK

  return ret_val;
  }


int write_image_colormap(const std::string & filename, const ImageGrid & data,
			 const Color::Colormap & colormap, 
			 const bool transpose_image,
			 const bool fixed_aspect_ratio,
			 const double horiz_vert_ratio)
  {
  int ret_val = -1;  // Error state

  if (filename.size() < 1)
    throw Except("Empty filename", __FILE__, __LINE__);

  const unsigned int num_pixels = data.get_num_pixels();

  if (num_pixels < 1)
    throw Except("No data for image", __FILE__, __LINE__);

  const std::vector<unsigned int> & image_dims = data.get_image_dims();

  // Determine which axis is which based on transpose_image flag

  const unsigned int image_horizontal_pixels 
    = transpose_image ? image_dims[1] : image_dims[0];

  const unsigned int image_vertical_pixels 
    = transpose_image ? image_dims[0] : image_dims[1];

  unsigned int horizontal_pixels = 0;
  unsigned int vertical_pixels = 0;

  unsigned int horizontal_offset = 0;
  unsigned int vertical_offset = 0;

  if (fixed_aspect_ratio)
    {
    if (!(horiz_vert_ratio > 0.0))
      throw Except("Non-positive aspect ratio", __FILE__, __LINE__);

    const double image_horiz_vert_ratio 
      = static_cast<double>(image_horizontal_pixels) 
        / static_cast<double>(image_vertical_pixels);

    if (image_horiz_vert_ratio > horiz_vert_ratio)
      {
      horizontal_pixels = image_horizontal_pixels;

      const double fl_vertical_pixels 
	= static_cast<double>(image_horizontal_pixels) / horiz_vert_ratio;

      vertical_pixels = static_cast<unsigned int>(fl_vertical_pixels);

      if (vertical_pixels < 1)
	vertical_pixels = 1;
      }
    else
      {
      vertical_pixels = image_vertical_pixels;

      const double fl_horizontal_pixels 
	= static_cast<double>(image_vertical_pixels) * horiz_vert_ratio;

      horizontal_pixels = static_cast<unsigned int>(fl_horizontal_pixels);

      if (horizontal_pixels < 1)
	horizontal_pixels = 1;
      }

    if (vertical_pixels > image_vertical_pixels)
      vertical_offset = (vertical_pixels - image_vertical_pixels) / 2;
    }
  else
    {
    horizontal_pixels = image_horizontal_pixels;
    vertical_pixels = image_vertical_pixels;
    }

  const Magick::Geometry geometry(horizontal_pixels, vertical_pixels);

  Magick::ColorRGB magick_color(0.0, 0.0, 0.0);

  Magick::Image image(geometry, magick_color);

  image.type(Magick::TrueColorMatteType);
  image.depth(8);

  Magick::Pixels pixel_cache(image);

  Magick::PixelPacket * image_pixels 
    = pixel_cache.get(horizontal_offset, vertical_offset, 
		      image_horizontal_pixels, image_vertical_pixels);

  // NOTE: Image is stored beginning with upper-left corner, proceding in rows

  std::vector<unsigned int> index(2);

  for (unsigned int j = image_vertical_pixels; j > 0; j--)
    {
    if (transpose_image)
      index[0] = j - 1;
    else
      index[1] = j - 1;

    for (unsigned int i = 0; i < image_horizontal_pixels; i++)
      {
      if (transpose_image)
	index[1] = i;
      else
	index[0] = i;

      const Color::Color & color = data[index];

      magick_color.red(color.get_red());
      magick_color.green(color.get_green());
      magick_color.blue(color.get_blue());

      // ImageMagick (prior to 6.3.0) defines alpha as level of transparency
      // ... Magick++ seems to follow this definition for now ...

      const double transparency = 1.0 - color.get_alpha();

      magick_color.alpha(transparency);

      *(image_pixels) = magick_color;

      image_pixels++;
      }
    }

  pixel_cache.sync();

  image.syncPixels();

  // Build the colormap legend

  const unsigned int bar_width = 32;

  const unsigned int bar_height 
    = std::min(vertical_pixels, 
	       std::max(static_cast<unsigned int>(4), 
			static_cast<unsigned int>((3 * vertical_pixels) / 4)));

  if ((bar_height > 20) && (horizontal_pixels >= (2 * bar_width)))
    {
    const unsigned int bar_edge = horizontal_pixels - (2 * bar_width);
    const unsigned int bar_space 
      = std::min((vertical_pixels / 8), ((vertical_pixels - bar_height) / 2));

    Magick::PixelPacket * bar_pixels 
      = pixel_cache.get(bar_edge, bar_space, bar_width, bar_height);

    const double min_value = colormap.get_min_value();
    const double max_value = colormap.get_max_value();

    const double value_range = max_value - min_value;

    const double value_delta 
      = value_range / static_cast<double>(bar_height - 1);

    for (unsigned int j = bar_height; j > 0; j--)
      {
      const double current_value 
	= min_value + (static_cast<double>(j) * value_delta);

      Color::Color current_color;

      colormap.get_color(current_value, current_color);

      magick_color.red(current_color.get_red());
      magick_color.green(current_color.get_green());
      magick_color.blue(current_color.get_blue());

      // ImageMagick (prior to 6.3.0) defines alpha as level of transparency
      // ... Magick++ seems to follow this definition for now ...

      const double transparency = 1.0 - current_color.get_alpha();

      magick_color.alpha(transparency);

      for (unsigned int i = 0; i < bar_width; i++)
	{
	(*bar_pixels) = magick_color;
	bar_pixels++;
	}
      }

    pixel_cache.sync();

    image.syncPixels();

    // Write some labels for the colormap legend

    image.penColor("white");

    image.fontPointsize(20);

    const unsigned int num_labels = 5;

    const double label_delta 
      = value_range / static_cast<double>(num_labels - 1);

    const unsigned int label_offset_delta 
      = std::max(static_cast<unsigned int>(bar_height / (num_labels - 1)), 
		 static_cast<unsigned int>(1));

    const unsigned int label_offset_start 
      = (bar_space > 10) ? (bar_space - 10) : 1;

    std::stringstream label_stream;

    for (unsigned int label_index = 0; label_index < num_labels; label_index++)
      {
      label_stream.str("");  // Clear the stream

      const double label_value 
	= max_value - (static_cast<double>(label_index) * label_delta);

      label_stream << label_value;

      const unsigned int label_start_horiz = bar_edge - 10;
      const unsigned int label_start_vert = 20;

      const unsigned int label_offset_horiz = 0;
      const unsigned int label_offset_vert 
	= label_offset_start  + (label_index * label_offset_delta);

      Magick::Geometry label_boundary(label_start_horiz, label_start_vert,
				      label_offset_horiz, label_offset_vert);

      image.annotate(label_stream.str().c_str(), label_boundary, 
		     Magick::EastGravity);
      }

    ret_val = 0;  // OK

    image.write(filename);
    }

  return ret_val;
  }


int write_colormap(const std::string & filename, 
		   const Color::Colormap & colormap, 
		   const unsigned int value_axis_length,
		   const unsigned int width, const bool horizontal,
		   const bool limit_values, 
		   const double min_legend_value, 
		   const double max_legend_value)
  {
  int ret_val = -1;  // Error state

  // Perform some checks

  if (value_axis_length < 2)
    throw Except("Improper colormap legend length", __FILE__, __LINE__);

  if (width < 1)
    throw Except("Improper colormap legend width", __FILE__, __LINE__);

  if (limit_values)
    if (min_legend_value > max_legend_value)
      throw Except("Incompatible value limits");  // Equal is OK

  // Build the ImageGrid object

  std::vector<unsigned int> grid_dims(2);

  grid_dims[0] = width;
  grid_dims[1] = value_axis_length;

  ImageGrid colormap_image(grid_dims);

  // Determine the range of colors
  
  const double min_value 
    = limit_values ? min_legend_value : colormap.get_min_value();

  const double max_value 
    = limit_values ? max_legend_value : colormap.get_max_value();

  const bool use_log_values = colormap.log_scale_active();

  const double proc_min_value 
    = use_log_values ? std::log10(min_value) : min_value;

  const double proc_max_value 
    = use_log_values ? std::log10(max_value) : max_value;

  const double delta_proc_value = (proc_max_value - proc_min_value) 
    / static_cast<double>(value_axis_length);
  
  // Start on a pixel center

  const double start_proc_value = proc_min_value + (0.5 * delta_proc_value);

  // Build the colormap image
  
  std::vector<unsigned int> pixel_index(2);

  for (unsigned int i = 0; i < width; i++)
    {
    double current_proc_value = start_proc_value;

    pixel_index[0] = i;

    for (unsigned int j = 0; j < value_axis_length; j++)
      {
      pixel_index[1] = j;

      const double current_value = use_log_values 
	? std::pow(10.0, current_proc_value) : current_proc_value;

      colormap.get_color(current_value, colormap_image[pixel_index]);

      current_proc_value += delta_proc_value;
      }
    }

  ret_val = write_image(filename, colormap_image, horizontal);  // HORIZ WORKS?

  return ret_val;
  }


int write_colormap(const std::string & colormap_image_filename, 
		   const std::string & colormap_info_filename,
		   const Color::Colormap & colormap, 
		   const unsigned int value_axis_length,
		   const unsigned int width, const bool horizontal,
		   const bool limit_values, 
		   const double min_legend_value, 
		   const double max_legend_value,
		   const bool use_comment_markers,
		   const char comment_marker)
  {
  int ret_val = write_colormap(colormap_image_filename, colormap, 
			       value_axis_length, width, horizontal, 
			       limit_values, min_legend_value, 
			       max_legend_value);

  // Use only the unqualified name for the colormap image

  std::string cmap_image_path;
  std::string cmap_image_name;

  Utils::split_name_path_file(colormap_image_filename,
			      cmap_image_path, cmap_image_name);

  if (ret_val >= 0)
    ret_val = write_colormap_info(cmap_image_name, 
				  colormap_info_filename, 
				  colormap, value_axis_length,
				  width, horizontal, limit_values, 
				  min_legend_value, max_legend_value,
				  use_comment_markers, comment_marker);

  return ret_val;
  }


int write_colormap_info(const std::string & colormap_image_filename,
			const std::string & colormap_info_filename,
			const Color::Colormap & colormap, 
			const unsigned int value_axis_length,
			const unsigned int width, const bool horizontal,
			const bool limit_values, 
			const double min_legend_value, 
			const double max_legend_value,
			const bool use_comment_markers,
			const char comment_marker)
  {
  using std::endl;

  int ret_val = -1;  // Error state

  std::ofstream f(colormap_info_filename.c_str(), std::ios::trunc);

  std::string line_prefix = "";

  if (use_comment_markers)
    line_prefix += comment_marker + " ";

  f << line_prefix << "ColormapImage " << colormap_image_filename << endl;
  f << line_prefix << "MinValue " << colormap.get_min_value() << endl;
  f << line_prefix << "MaxValue " << colormap.get_max_value() << endl;
  f << line_prefix << "ValueAxisLength " << value_axis_length << endl;
  f << line_prefix << "WidthPixels " << width << endl;

  f << line_prefix << "HorizontalColormap ";

  if (horizontal)
    f << "true";
  else
    f << "false";

  f << endl;

  f << line_prefix << "LimitLegendValues ";

  if (limit_values)
    f << "true";
  else
    f << "false";
  
  f << endl;

  if (limit_values)
    {
    f << "MinLegendValue " << min_legend_value << endl;
    f << "MaxLegendValue " << max_legend_value << endl;
    }

  f << line_prefix << "UseLogValues ";

  if (colormap.log_scale_active())
    f << "true";
  else
    f << "false";

  f << endl;

  f << line_prefix << "ColormapData" << endl;

  const unsigned int num_entries = colormap.get_num_entries();

  for (unsigned int index = 0; index < num_entries; index++)
    {
    const double current_value = colormap.get_index_value(index);
    const Color::Color & current_color = colormap[index];

    const double red = current_color.get_red();
    const double green = current_color.get_green();
    const double blue = current_color.get_blue();

    const double alpha = current_color.get_alpha();

    f << current_value << " " << red << " " << green << " " << blue
      << " " << alpha << endl;
    }

  ret_val = 0;  // OK, I guess ...

  return ret_val;
  }


int write_colormap_info(const std::string & colormap_info_filename,
			const Color::Colormap & colormap, 
			const unsigned int value_axis_length,
			const unsigned int width, const bool horizontal,
			const bool limit_values, 
			const double min_legend_value, 
			const double max_legend_value,
			const bool use_comment_markers,
			const char comment_marker)
  {
  using std::endl;

  int ret_val = -1;  // Error state

  std::ofstream f(colormap_info_filename.c_str(), std::ios::trunc);

  std::string line_prefix = "";

  if (use_comment_markers)
    line_prefix += comment_marker + " ";

  f << line_prefix << "MinValue " << colormap.get_min_value() << endl;
  f << line_prefix << "MaxValue " << colormap.get_max_value() << endl;
  f << line_prefix << "ValueAxisLength " << value_axis_length << endl;
  f << line_prefix << "WidthPixels " << width << endl;

  f << line_prefix << "HorizontalColormap ";

  if (horizontal)
    f << "true";
  else
    f << "false";

  f << endl;

  f << line_prefix << "LimitLegendValues ";

  if (limit_values)
    f << "true";
  else
    f << "false";
  
  f << endl;

  if (limit_values)
    {
    f << "MinLegendValue " << min_legend_value << endl;
    f << "MaxLegendValue " << max_legend_value << endl;
    }

  f << line_prefix << "UseLogValues ";

  if (colormap.log_scale_active())
    f << "true";
  else
    f << "false";

  f << endl;

  f << line_prefix << "ColormapData" << endl;

  const unsigned int num_entries = colormap.get_num_entries();

  for (unsigned int index = 0; index < num_entries; index++)
    {
    const double current_value = colormap.get_index_value(index);
    const Color::Color & current_color = colormap[index];

    const double red = current_color.get_red();
    const double green = current_color.get_green();
    const double blue = current_color.get_blue();

    const double alpha = current_color.get_alpha();

    f << current_value << " " << red << " " << green << " " << blue
      << " " << alpha << endl;
    }

  ret_val = 0;  // OK, I guess ...

  return ret_val;
  }


}  // End namespace ImageMagick
}  // End namespace Image
}  // End namespace QuickFlash


#endif  // USE_MAGICK
