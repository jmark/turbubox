// C++ header file quickflash_color_color.hpp

/*
  By Nathan C. Hearn
    February 11, 2008

  Color-related classes.
*/


#ifndef QUICKFLASH_COLOR_COLOR_HPP
#define QUICKFLASH_COLOR_COLOR_HPP


namespace QuickFlash
{
namespace Color
{

enum ColorChannel { Red, Green, Blue, Alpha };


class Color
  {
  public :

    Color() : 
      red_value(0.0), green_value(0.0), blue_value(0.0), alpha_value(0.0)
      { }

    Color(const double gray, const double alpha=1.0) :
      red_value(0.0), green_value(0.0), blue_value(0.0), alpha_value(0.0)
      { reset(gray, alpha); }


    Color(const double red, const double green, const double blue,
	  const double alpha=1.0) :
      red_value(0.0), green_value(0.0), blue_value(0.0), alpha_value(0.0)
      { reset(red, green, blue, alpha); }


    Color(const Color & source) :
      red_value(0.0), green_value(0.0), blue_value(0.0), alpha_value(0.0)
      { reset(source); }


    ~Color() { }

    Color & operator=(const Color & source)
      {
      reset(source);
      return *this;
      }

    Color & operator+=(const Color & source)
      {
      blend_color(source);
      return *this;
      }

    void reset();

    void reset(const double gray, const double alpha=1.0)
      { reset(gray, gray, gray, alpha); }

    void reset(const double red, const double green,
	       const double blue, const double alpha=1.0);

    void reset(const Color & source);

    void reset(const Color & color1, const Color & color2, 
	       const double color1_fract=0.5);

    double get_red() const { return red_value; }
    double get_green() const { return green_value; }
    double get_blue() const { return blue_value; }

    double get_alpha() const { return alpha_value; }

    double get_channel(const ColorChannel channel) const;

    void set_gray(const double value)
      {
      set_red(value);
      set_green(value);
      set_blue(value);
      }

    void set_red(const double value);
    void set_green(const double value);
    void set_blue(const double value);

    void set_alpha(const double value);

    void set_channel(const ColorChannel channel, const double value);

    void set_color_hsv(const double hue, const double saturation,
		       const double value, const double alpha=1.0);

    void blend_color(const Color & source, const double local_scale=1.0,
		     const double source_scale=1.0);

    void overlay_color(const Color & source, const bool add_below=false);

  private :
    
    void limit_values();

  private :

    double red_value;
    double green_value;
    double blue_value;

    double alpha_value;  // Level of opacity (1.0 = fully opaque)
  };


class IntColor
  {
  public :

    IntColor() : 
      red_value(0), green_value(0), blue_value(0), alpha_value(0)
      { }

    IntColor(const unsigned int red, const unsigned int green,
	     const unsigned int blue, const unsigned int alpha) :
      red_value(0), green_value(0), blue_value(0), alpha_value(0)
      { reset(red, green, blue, alpha); }

    IntColor(const IntColor & source) :
      red_value(0), green_value(0), blue_value(0), alpha_value(0)
      { reset(source); }

    ~IntColor() { }

    IntColor & operator=(const IntColor & source)
      {
      reset(source);
      return *this;
      }

    void reset()
      {
      red_value = 0;
      green_value = 0;
      blue_value = 0;

      alpha_value = 0;
      }

    void reset(const unsigned int red, const unsigned int green,
	       const unsigned int blue, const unsigned int alpha)
      {
      red_value = red;
      green_value = green;
      blue_value = blue;
      
      alpha_value = alpha;
      }

    void reset(const IntColor & source)
      {
      if (&source !=this)
	{
	red_value = source.red_value;
	green_value = source.green_value;
	blue_value = source.blue_value;

	alpha_value = source.alpha_value;
	}
      }

    unsigned int get_red() const { return red_value; }
    unsigned int get_green() const { return green_value; }
    unsigned int get_blue() const { return blue_value; }

    unsigned int get_alpha() const { return alpha_value; }

    void set_red(const unsigned int value) { red_value = value; }
    void set_green(const unsigned int value) { green_value = value; }
    void set_blue(const unsigned int value) { blue_value = value; }

    void set_alpha(const unsigned int value) { alpha_value = value; }

  private :

    unsigned int red_value;
    unsigned int green_value;
    unsigned int blue_value;

    unsigned int alpha_value;
  };


class IntColorConverter
  {
  public :
    
    IntColorConverter(const unsigned int levels=256) :
      max_level(0), oneOver_levels(0.0)
      { reset(levels); }

    IntColorConverter(const IntColorConverter & source) :
      max_level(0), oneOver_levels(0.0)
      { reset(source); }

    ~IntColorConverter() { }

    IntColorConverter & operator=(const IntColorConverter & source)
      {
      reset(source);
      return *this;
      }

    void reset(const unsigned int levels=256)
      {
      max_level = (levels > 0) ? (levels - 1) : 0;

      oneOver_levels = 1.0 / static_cast<double>(max_level + 1);
      }

    void reset(const IntColorConverter & source)
      {
      if (&source != this)
	{
	max_level = source.max_level;
	oneOver_levels = source.oneOver_levels;
	}
      }

    void get_int_color(const Color & color, IntColor & int_color) const;

  private :

    unsigned int max_level;

    double oneOver_levels;
  };


}  // End namespace Color
}  // End namespace QuickFlash


#endif  // QUICKFLASH_COLOR_COLOR_HPP
