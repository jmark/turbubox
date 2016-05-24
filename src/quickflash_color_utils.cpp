// C++ program file quickflash_color_utils.cpp

/*
  By Nathan C. Hearn
     June 14, 2008

  Color utility functions.
*/


#include "quickflash_color_utils.hpp"


#include <cmath>
#include <algorithm>


namespace QuickFlash
{
namespace Color
{

void convert_offset_hue_channel(const double offset_hue, double & channel)
  {
  const double one_sixth = 1.0 / 6.0;
  const double one_third = 1.0 / 3.0;

  // Here, the center of the channel peak is at offset_hue == 0

  double work_hue = offset_hue;

  // Rotate the offset hue so that it is in the range [-0.5,0.5]

  work_hue += 0.5;  // Work in range [0,1]

  if ((work_hue < 0.0) || (work_hue > 1.0))
    {
    work_hue -= std::floor(work_hue);

    if (work_hue < 0.0)
      work_hue = 0.0;
    else if (work_hue > 1.0)
      work_hue = 1.0;
    }

  work_hue -= 0.5;  // Back to range [-0.5, 0.5]

  // Determine the channel level

  const double hue_dist = std::fabs(work_hue);

  // Here, we set the channel value as follows:
  //
  //   hue_dist < 1/6 -> value = 1
  //   hue_dist > 1/3 -> value = 0
  //   ... otherwise  -> value = 1 - (6 * (hue_dist - 1/6))

  double work_channel = 0.0;  // Default

  if (hue_dist < one_sixth)
    work_channel = 1.0;
  else if (hue_dist < one_third)
    work_channel = 1.0 - (6.0 * (hue_dist - one_sixth));

  if (work_channel < 0.0)
    work_channel = 0.0;
  else if (work_channel > 1.0)
    work_channel = 1.0;

  channel = work_channel;
  }


void convert_hue_rgb(const double hue, double & red, double & green, 
		     double & blue)
  {
  const double one_third = 1.0 / 3.0;
  const double two_thirds = 2.0 / 3.0;

  // Determine normalized channel values for the given hue

  double work_red = 0.0;
  double work_green = 0.0;
  double work_blue = 0.0;

  // Red hue is centered at hue = 0
  // Green hue is centered at hue = 1/3
  // Blue hue is centered at hue = 2/3

  // Center the hues

  convert_offset_hue_channel(hue, work_red);
  convert_offset_hue_channel((hue - one_third), work_green);
  convert_offset_hue_channel((hue - two_thirds), work_blue);

  // Find the maximum

  const double max_value = std::max(work_red, std::max(work_green, work_blue));

  // Scale all colors so that the max value is one

  const double channel_scale = 1.0 / max_value;  // ERROR CHECKING???

  red = work_red * channel_scale;
  green = work_green * channel_scale;
  blue = work_blue * channel_scale;
  }


void convert_hsv_rgb(const double hue, const double saturation, 
		     const double value, double & red, double & green, 
		     double & blue)
  {
  // ERROR CHECKING???

  double work_red = 0.0;
  double work_green = 0.0;
  double work_blue = 0.0;

  convert_hue_rgb(hue, work_red, work_green, work_blue);

  if (saturation < 1.0)
    {
    // ERROR CHECKING???

    const double sat_scale = 1.0 - saturation;
	
    work_red += sat_scale * (1.0 - work_red);
    work_green += sat_scale * (1.0 - work_green);
    work_blue += sat_scale * (1.0 - work_blue);
    }

  // Scale by the value

  // ERROR CHECKING???

  work_red *= value;
  work_green *= value;
  work_blue *= value;

  if (work_red < 0.0)
    work_red = 0.0;
  else if (work_red > 1.0)
    work_red = 1.0;

  if (work_green < 0.0)
    work_green = 0.0;
  else if (work_green > 1.0)
    work_green = 1.0;

  if (work_blue < 0.0)
    work_blue = 0.0;
  else if (work_blue > 1.0)
    work_blue = 1.0;

  red = work_red;
  green = work_green;
  blue = work_blue;
  }


}  // End namespace Color
}  // End namespace QuickFlash
