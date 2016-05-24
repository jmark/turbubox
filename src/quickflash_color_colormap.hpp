// C++ header file quickflash_color_colormap.hpp

/*
  By Nathan C. Hearn
    February 11, 2008

  Colormap-related classes.
*/


#ifndef QUICKFLASH_COLOR_COLORMAP_HPP
#define QUICKFLASH_COLOR_COLORMAP_HPP


#include <string>
#include <list>
#include <vector>
#include "quickflash_color_color.hpp"


namespace QuickFlash
{
namespace Color
{

// Class ValueColor

class ValueColor : public Color
  {
  public :

    ValueColor() : Color(), map_value(0.0) { }

    ValueColor(const double value, const double gray, const double alpha=1.0) :
      Color(), map_value(0.0)
      { reset(value, gray, alpha); }

    ValueColor(const double value, const double red, 
	       const double green, const double blue,
	       const double alpha=1.0) :
      Color(), map_value(0.0)
      { reset(value, red, green, blue, alpha); }

    ValueColor(const double value, const Color & color) :
      Color(), map_value(0.0)
      { reset(value, color); }

    ValueColor(const ValueColor & source) :
      Color(), map_value(0.0)
      { reset(source); }

    ~ValueColor() { }

    ValueColor & operator=(const ValueColor & source)
      {
      reset(source);
      return *this;
      }

    ValueColor & operator=(const Color & source)
      {
      reset_color(source);  // NOTE: We don't change the value component
      return *this;
      }

    bool operator<(const ValueColor & comp_color) const
      { return (map_value < comp_color.map_value); }

    void reset()
      {
      map_value = 0.0;
      Color::reset();
      }

    void reset(const double value, const double gray, const double alpha=1.0)
      {
      map_value = value;
      Color::reset(gray, alpha);
      }

    void reset(const double value, const double red, 
	       const double green, const double blue,
	       const double alpha=1.0)
      {
      map_value = value;
      Color::reset(red, green, blue, alpha);
      }

    void reset(const double value, const Color & color)
      {
      map_value = value;
      reset_color(color);
      }

    void reset(const ValueColor & source)
      {
      if (&source != this)
	{
	map_value = source.map_value;
	reset_color(source);
	}
      }

    void reset_color(const Color & color)
      { Color::reset(color); }

    double get_value() const { return map_value; }

    void set_value(const double value) { map_value = value; }

  private :

    double map_value;
  };


// Class Colormap

class Colormap
  {
  public :

    Colormap();

    Colormap(const double min_value, const double max_value, 
	     const bool log_scale=false, const double opacity=1.0);

    Colormap(const std::list<double> & values, 
	     const std::list<Color> & colors, const bool log_scale=false);

    Colormap(const std::list<ValueColor> & value_colors, 
	     const bool log_scale=false);

    Colormap(const Colormap & source);

    ~Colormap() { }

    Colormap & operator=(const Colormap & source)
      {
      reset(source);
      return *this;
      }

    const Color & operator[](const unsigned int index) const
      { return get_index_color(index); }

    void reset();

    void reset(const double min_value, const double max_value, 
	       const bool log_scale=false, const double opacity=1.0);

    void reset(const std::list<double> & values, 
	       const std::list<Color> & colors, const bool log_scale=false);

    void reset(const std::list<ValueColor> & value_colors, 
	       const bool log_scale=false);

    void reset(const Colormap & source);

    void rescale_values(const double new_min_value, const double new_max_value,
			const bool log_scale=false);

    double get_min_value() const { return min_value; }
    double get_max_value() const { return max_value; }

    bool map_initialized() const 
      { 
      bool ret_val = false;

      if (num_entries > 1)
	ret_val = true;

      return ret_val;
      }

    bool log_scale_active() const { return use_log_values; }

    unsigned int get_num_entries() const { return num_entries; }

    void get_color(const double value, Color & color) const;

    double get_index_value(const unsigned int index) const
      { return value_list[index]; }

    const Color & get_index_color(const unsigned int index) const
      { return color_list[index]; }

    const std::vector<double> & get_value_list() const { return value_list; }

    const std::vector<Color> & get_color_list() const { return color_list; }

    void write_colormap_info(const std::string & filename,
			     const bool use_comment_markers=true,
			     const char comment_marker='#') const;
    
  private :

    void set_colors(const std::list<ValueColor> & sorted_value_colors,
		    const bool log_scale=false);

    void locate_value(const double processed_value, 
		      unsigned int & before_index, unsigned int & after_index,
		      double & value_before_delta, double & interval_size) 
      const;

  private :

    bool use_log_values;

    unsigned int num_entries;

    std::vector<double> value_list;
    std::vector<double> processed_value_list;

    double min_value;
    double max_value;

    std::vector<Color> color_list;
  };


// Utility functions

int extract_numbers(const std::string & line, 
		    std::vector<double> & number_list);

int load_colormap(const std::string & filename, Colormap & colormap,
		  const bool log_scale=false);


}  // End namespace Color
}  // End namespace QuickFlash


#endif  // QUICKFLASH_COLOR_COLORMAP_HPP
