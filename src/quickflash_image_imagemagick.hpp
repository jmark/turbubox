// C++ header file quickflash_image_imagemagick.hpp

/*
  By Nathan C. Hearn
     February 12, 2008

  Routines for converting ImageGrid object data into bitmapped images.

  Requires ImageMagick's Magick++ bindings.
*/


#ifndef QUICKFLASH_IMAGE_IMAGEMAGICK_HPP
#define QUICKFLASH_IMAGE_IMAGEMAGICK_HPP

#ifdef USE_MAGICK


#include <string>
#include "quickflash_image_imagegrid.hpp"
#include "quickflash_color_colormap.hpp"


namespace QuickFlash
{
namespace Image
{
namespace ImageMagick
{

int write_image(const std::string & filename, const ImageGrid & data,
		const bool transpose_image=false);

int write_image_colormap(const std::string & filename, const ImageGrid & data,
			 const Color::Colormap & colormap, 
			 const bool transpose_image=false,
			 const bool fixed_aspect_ratio=false,
			 const double horiz_vert_ratio=(4.0/3.0));

int write_colormap(const std::string & filename, 
		   const Color::Colormap & colormap, 
		   const unsigned int value_axis_length,
		   const unsigned int width, const bool horizontal=false,
		   const bool limit_values=false, 
		   const double min_legend_value=0.0, 
		   const double max_legend_value=1.0);

int write_colormap(const std::string & colormap_image_filename, 
		   const std::string & colormap_info_filename,
		   const Color::Colormap & colormap, 
		   const unsigned int value_axis_length,
		   const unsigned int width, const bool horizontal=false,
		   const bool limit_values=false, 
		   const double min_legend_value=0.0, 
		   const double max_legend_value=1.0,
		   const bool use_comment_markers=false,
		   const char comment_marker='#');

int write_colormap_info(const std::string & colormap_image_filename, 
			const std::string & colormap_info_filename,
			const Color::Colormap & colormap, 
			const unsigned int value_axis_length,
			const unsigned int width, const bool horizontal=false,
			const bool limit_values=false, 
			const double min_legend_value=0.0, 
			const double max_legend_value=1.0,
			const bool use_comment_markers=false,
			const char comment_marker='#');

int write_colormap_info(const std::string & colormap_info_filename,
			const Color::Colormap & colormap, 
			const unsigned int value_axis_length,
			const unsigned int width, const bool horizontal=false,
			const bool limit_values=false, 
			const double min_legend_value=0.0, 
			const double max_legend_value=1.0,
			const bool use_comment_markers=false,
			const char comment_marker='#');


}  // End namespace ImageMagick
}  // End namespace Image
}  // End namespace QuickFlash


#endif  // USE_MAGICK


#endif  // QUICKFLASH_IMAGE_IMAGEMAGICK_HPP
