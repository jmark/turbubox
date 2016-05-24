// C++ program file quickflash_color_color.cpp

/*
  By Nathan C. Hearn
    February 11, 2008

  Color-related classes.
*/


#include "quickflash_color_color.hpp"
#include "quickflash_color_utils.hpp"


namespace QuickFlash
{
namespace Color
{

// Class Color

void Color::reset()
  {
  red_value = 0.0;
  green_value = 0.0;
  blue_value = 0.0;

  alpha_value = 0.0;
  }


void Color::reset(const double red, const double green, const double blue, 
		  const double alpha)
  {
  red_value = red;
  green_value = green;
  blue_value = blue;

  alpha_value = alpha;

  limit_values();
  }


void Color::reset(const Color & source)
  {
  if (&source != this)
    {
    red_value = source.red_value;
    green_value = source.green_value;
    blue_value = source.blue_value;

    alpha_value = source.alpha_value;
    }
  }


void Color::reset(const Color & color1, const Color & color2, 
		  const double color1_fract)
  {
  reset(color1);

  const double color2_fract = 1.0 - color1_fract;  // ERROR CHECKING???

  blend_color(color2, color1_fract, color2_fract);
  }


double Color::get_channel(const ColorChannel channel) const
  {
  double ret_val = 0.0;

  switch (channel)
    {
    case Red :
      ret_val = red_value;
      break;

    case Green :
      ret_val = green_value;
      break;
	  
    case Blue :
      ret_val = blue_value;
      break;

    case Alpha :
      ret_val = alpha_value;
      break;
    }

  return ret_val;
  }


void Color::set_red(const double value)
  {
  red_value = value;

  if (red_value < 0.0)
    red_value = 0.0;
  else if (red_value > 1.0)
    red_value = 1.0;
  }


void Color::set_green(const double value)
  {
  green_value = value;

  if (green_value < 0.0)
    green_value = 0.0;
  else if (green_value > 1.0)
    green_value = 1.0;
  }


void Color::set_blue(const double value)
  {
  blue_value = value;

  if (blue_value < 0.0)
    blue_value = 0.0;
  else if (blue_value > 1.0)
    blue_value = 1.0;
  }


void Color::set_alpha(const double value)
  {
  alpha_value = value;

  if (alpha_value < 0.0)
    alpha_value = 0.0;
  else if (alpha_value > 1.0)
    alpha_value = 1.0;
  }


void Color::set_channel(const ColorChannel channel, const double value)
  {
  double norm_value = value;

  if (norm_value < 0.0)
    norm_value = 0.0;
  else if (norm_value > 1.0)
    norm_value = 1.0;

  switch (channel)
    {
    case Red :
      red_value = norm_value;
      break;

    case Green :
      green_value = norm_value;
      break;
	  
    case Blue :
      blue_value = norm_value;
      break;

    case Alpha :
      alpha_value = norm_value;
      break;
    }
  }


void Color::set_color_hsv(const double hue, const double saturation,
			  const double value, const double alpha)
  {
  convert_hsv_rgb(hue, saturation, value, red_value, green_value, blue_value);

  alpha_value = alpha;
  }


void Color::blend_color(const Color & source, const double local_scale,
			const double source_scale)
  {
  // Just add the color according to the scales

  red_value *= local_scale;
  red_value += source.red_value * source_scale;

  green_value *= local_scale;
  green_value += source.green_value * source_scale;

  blue_value *= local_scale;
  blue_value += source.blue_value * source_scale;

  // I'm less sure what to do with alpha ... just treat the same for now

  alpha_value *= local_scale;
  alpha_value += source.alpha_value * source_scale;

  limit_values();
  }


void Color::overlay_color(const Color & source, const bool add_below)
  {
  // Determine the amount of the resulting color that is due to the this 
  // color (local) and that due to the source color.  If the source is added
  // below this color, use this color's opacity (alpha) for the local 
  // fraction; otherwise use the source's transparency (complement of alpha) 
  // for the local fraction.

  const double local_color_fract 
    = add_below ? alpha_value : (1.0 - source.alpha_value);

  const double source_color_fract = 1.0 - local_color_fract;

  red_value *= local_color_fract;
  red_value += source.red_value * source_color_fract;

  green_value *= local_color_fract;
  green_value += source.green_value * source_color_fract;

  blue_value *= local_color_fract;
  blue_value += source.blue_value * source_color_fract;

  // To reflect the increase in opacity, take the product of the transparencies

  const double local_trans = 1.0 - alpha_value;          // 1.0 = transparent
  const double source_trans = 1.0 - source.alpha_value;  // ...

  const double combined_trans = local_trans * source_trans;

  alpha_value = 1.0 - combined_trans;  // Convert back to opacity

  // Just in case ... limit the values

  limit_values();
  }


void Color::limit_values()
  {
  if (red_value < 0.0)
    red_value = 0.0;
  else if (red_value > 1.0)
    red_value = 1.0;

  if (green_value < 0.0)
    green_value = 0.0;
  else if (green_value > 1.0)
    green_value = 1.0;

  if (blue_value < 0.0)
    blue_value = 0.0;
  else if (blue_value > 1.0)
    blue_value = 1.0;

  if (alpha_value < 0.0)
    alpha_value = 0.0;
  else if (alpha_value > 1.0)
    alpha_value = 1.0;
  }


// Class IntColorConverter

void IntColorConverter::get_int_color(const Color & color, 
				      IntColor & int_color) const
  {
  unsigned int red_value 
    = static_cast<unsigned int>(color.get_red() * oneOver_levels);

  if (red_value > max_level)
    red_value = max_level;

  unsigned int green_value 
    = static_cast<unsigned int>(color.get_green() * oneOver_levels);

  if (green_value > max_level)
    green_value = max_level;

  unsigned int blue_value 
    = static_cast<unsigned int>(color.get_blue() * oneOver_levels);

  if (blue_value > max_level)
    blue_value = max_level;

  unsigned int alpha_value 
    = static_cast<unsigned int>(color.get_alpha() * oneOver_levels);

  if (alpha_value > max_level)
    alpha_value = max_level;

  int_color.reset(red_value, green_value, blue_value, alpha_value);
  }


}  // End namespace Color
}  // End namespace QuickFlash
