// C++ header file quickflash_color_utils.hpp

/*
  By Nathan C. Hearn
     June 14, 2008

  Color utility functions.
*/


#ifndef QUICKFLASH_COLOR_UTILS_HPP
#define QUICKFLASH_COLOR_UTILS_HPP


namespace QuickFlash
{
namespace Color
{

void convert_offset_hue_channel(const double offset_hue, double & channel);


void convert_hue_rgb(const double hue, double & red, double & green, 
		     double & blue);


void convert_hsv_rgb(const double hue, const double saturation, 
		     const double value, double & red, double & green, 
		     double & blue);


}  // End namespace Color
}  // End namespace QuickFlash


#endif  // QUICKFLASH_COLOR_UTILS_HPP
