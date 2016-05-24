// C++ program file quickflash_color_colormap.cpp

/*
  By Nathan C. Hearn
    February 11, 2008

  Colormap-related classes.
*/


#include "quickflash_color_colormap.hpp"
#include <cmath>
#include <list>
#include <vector>
#include <fstream>
#include <string>
#include "quickflash_except.hpp"
#include "quickflash_color_color.hpp"
#include "quickflash_utils_text.hpp"


namespace QuickFlash
{
namespace Color
{
// Class Colormap

Colormap::Colormap() : 
  use_log_values(false), num_entries(0), value_list(), processed_value_list(), 
  min_value(0.0), max_value(0.0), color_list() 
  { }


Colormap::Colormap(const double min_value, const double max_value, 
		   const bool log_scale, const double opacity) : 
  use_log_values(false), num_entries(0), value_list(), processed_value_list(), 
  min_value(0.0), max_value(0.0), color_list() 
  { reset(min_value, max_value, log_scale, opacity); }


Colormap::Colormap(const std::list<double> & values, 
		   const std::list<Color> & colors, const bool log_scale) : 
  use_log_values(false), num_entries(0), value_list(), processed_value_list(), 
  min_value(0.0), max_value(0.0), color_list() 
  { reset(values, colors, log_scale); }


Colormap::Colormap(const std::list<ValueColor> & value_colors, 
		   const bool log_scale) : 
  use_log_values(false), num_entries(0), value_list(), processed_value_list(), 
  min_value(0.0), max_value(0.0), color_list() 
  { reset(value_colors, log_scale); }



Colormap::Colormap(const Colormap & source) : 
  use_log_values(false), num_entries(0), value_list(), processed_value_list(), 
  min_value(0.0), max_value(0.0), color_list() 
  { reset(source); }


void Colormap::reset()
  {
  use_log_values = false;

  num_entries = 0;

  value_list.clear();
  processed_value_list.clear();

  min_value = 0.0;
  max_value = 0.0;

  color_list.clear();
  }


void Colormap::reset(const double min_value, const double max_value, 
		     const bool log_scale, const double opacity)
  {
  // ERROR CHECKING???

  const ValueColor start_color(min_value, 0.0, opacity);
  const ValueColor end_color(max_value, 1.0, opacity);

  std::list<ValueColor> value_color_list;

  value_color_list.push_back(start_color);
  value_color_list.push_back(end_color);

  value_color_list.sort();

  set_colors(value_color_list, log_scale);
  }


void Colormap::reset(const std::list<double> & values, 
		     const std::list<Color> & colors, 
		     const bool log_scale)
  {
  const unsigned int list_size = values.size();

  if (colors.size() != list_size)
    throw Except("Incompatible list sizes", __FILE__, __LINE__);

  std::list<ValueColor> value_color_list;

  std::list<double>::const_iterator value_iter = values.begin();
  const std::list<double>::const_iterator end_value_iter = values.end();

  std::list<Color>::const_iterator color_iter = colors.begin();

  while (value_iter != end_value_iter)
    {
    const ValueColor value_color(*value_iter, *color_iter);

    value_color_list.push_back(value_color);

    ++value_iter;
    ++color_iter;
    }

  value_color_list.sort();

  set_colors(value_color_list, log_scale);
  }


void Colormap::reset(const std::list<ValueColor> & value_colors, 
		     const bool log_scale)
  {
  std::list<ValueColor> value_color_list = value_colors;

  value_color_list.sort();

  set_colors(value_color_list, log_scale);
  }


void Colormap::reset(const Colormap & source)
  {
  if (&source != this)
    {
    use_log_values = source.use_log_values;

    num_entries = source.num_entries;

    value_list = source.value_list;
    processed_value_list = source.processed_value_list;

    min_value = source.min_value;
    max_value = source.max_value;

    color_list = source.color_list;
    }
  }


void Colormap::rescale_values(const double new_min_value, 
			      const double new_max_value,
			      const bool log_scale)
  {
  if (!map_initialized())
    throw Except("Colormap not initialized", __FILE__, __LINE__);

  // Check the entries

  if (!(new_min_value < new_max_value))
    throw Except("Min value must be less than max value", __FILE__, __LINE__);

  if (log_scale)
    if (!(new_min_value > 0.0))
      throw Except("Non-positive values not allowed for log scale", __FILE__, 
		   __LINE__);

  // Rescale the old values

  const double scale_factor 
    = (new_max_value - new_min_value) / (max_value - min_value);

  for (unsigned int index = 0; index < num_entries; index++)
    {
    const double old_value = value_list[index];

    const double new_value 
      = ((old_value - min_value) * scale_factor) + new_min_value;

    value_list[index] = new_value;

    processed_value_list[index] 
      = log_scale ? std::log10(new_value) : new_value;
    }

  // Reset the min and max, as well as the log_values boolean

  min_value = new_min_value;
  max_value = new_max_value;

  use_log_values = log_scale;
  }


void Colormap::get_color(const double value, Color & color) const
  {
  if (num_entries < 2)
    throw Except("Colormap not initialized", __FILE__, __LINE__);

  if (!(value > min_value))
    color = color_list[0];
  else if (!(value < max_value))
    color = color_list[num_entries - 1];
  else
    {
    const double processed_value = use_log_values ? std::log10(value) : value;

    unsigned int before_index = 0;
    unsigned int after_index = 0;

    double value_before_delta = 0.0;
    double interval_size = 0.0;

    locate_value(processed_value, before_index, after_index, 
		 value_before_delta, interval_size);

    if (!(interval_size > 0.0))
      throw Except("Zero interval size found", __FILE__, __LINE__);

    // The fraction of the color at before_index is given by
    // [ 1 - (value_before_delta / interval_size) ].

    double before_fract = 1.0 - (value_before_delta / interval_size);

    if (before_fract < 0.0)
      before_fract = 0.0;
    else if (before_fract > 1.0)
      before_fract = 1.0;

    color.reset(color_list[before_index], color_list[after_index],
		before_fract);
    }
  }


void Colormap::write_colormap_info(const std::string & filename,
				   const bool use_comment_markers,
				   const char comment_marker) const
  {
  using std::endl;

  if (num_entries < 2)
    throw Except("Colormap not initialized", __FILE__, __LINE__);

  std::ofstream f(filename.c_str(), std::ios::trunc);

  std::string line_prefix = "";

  if (use_comment_markers)
    line_prefix += comment_marker + " ";

  f << line_prefix << "Colormap " << filename << endl;
  f << line_prefix << "MinValue " << min_value << endl;
  f << line_prefix << "MaxValue " << max_value << endl;

  f << line_prefix << "UseLogValues ";

  if (use_log_values)
    f << "true";
  else
    f << "false";

  f << endl;

  f << line_prefix << "Data" << endl;

  for (unsigned int index = 0; index < num_entries; index++)
    {
    const Color & current_color = color_list[index];

    const double red = current_color.get_red();
    const double green = current_color.get_green();
    const double blue = current_color.get_blue();

    const double alpha = current_color.get_alpha();

    f << value_list[index] << " " << red << " " << green << " " << blue
      << " " << alpha << endl;
    }
  }


void Colormap::set_colors(const std::list<ValueColor> & sorted_value_colors,
			  const bool log_scale)
  {
  const unsigned int num_colors = sorted_value_colors.size();

  if (num_colors < 2)
    throw Except("Colormap requires at least two colors", __FILE__, __LINE__);

  value_list.resize(num_colors);
  color_list.resize(num_colors);

  std::list<ValueColor>::const_iterator color_iter 
    = sorted_value_colors.begin();

  const std::list<ValueColor>::const_iterator end_color_iter 
    = sorted_value_colors.end();

  unsigned int index = 0;

  while (color_iter != end_color_iter)
    {
    const ValueColor & current_color = *color_iter;

    value_list[index] = current_color.get_value();
    color_list[index] = current_color;
 
    index++;
    ++color_iter;
    }

  // Scan the values

  double last_value = value_list[0];

  for (unsigned int index = 1; index < num_colors; index++)
    {
    const double current_value = value_list[index];

    if (!(current_value > last_value))
      throw Except("Values out of order or duplicated", __FILE__, __LINE__);

    last_value = current_value;
    }

  num_entries = num_colors;

  min_value = value_list[0];
  max_value = value_list[num_colors - 1];

  // Compute logs if desired

  use_log_values = log_scale;

  if (log_scale)
    {
    if (!(min_value > 0.0))
      throw Except("Log scale requires positive values", __FILE__, __LINE__);

    processed_value_list.resize(num_colors);

    for (unsigned int index = 0; index < num_colors; index++)
      {
      const double log_value = std::log10(value_list[index]);

      processed_value_list[index] = log_value;
      }

    // Now that logs were used, double-check for duplicate values

    double last_processed_value = processed_value_list[0];

    for (unsigned int index = 1; index < num_colors; index++)
      {
      const double current_value = processed_value_list[index];

      if (!(current_value > last_processed_value))
	throw Except("Values out of order or duplicated", __FILE__, __LINE__);

      last_processed_value = current_value;
      }
    }
  else
    processed_value_list = value_list;
  }


void Colormap::locate_value(const double processed_value, 
			    unsigned int & before_index, 
			    unsigned int & after_index,
			    double & value_before_delta, 
			    double & interval_size) const
  {
  // Bound the value

  unsigned int lower_index = 0;
  unsigned int upper_index = num_entries - 1;

  while (lower_index < (upper_index - 1))
    {
    const unsigned int mid_index = (lower_index + upper_index) / 2;

    const double mid_value = processed_value_list[mid_index];

    if (mid_value > processed_value)
      upper_index = mid_index;
    else
      lower_index = mid_index;
    }

  before_index = lower_index;
  after_index = upper_index;

  const double lower_value = processed_value_list[lower_index];
  const double upper_value = processed_value_list[upper_index];

  value_before_delta = processed_value - lower_value;

  interval_size = upper_value - lower_value;
  }


// Utility functions

int extract_numbers(const std::string & line, 
		    std::vector<double> & number_list)
  {
  int ret_val = -1;  // Error state

  number_list.clear();

  std::string::size_type comment_pos = line.find('#');

  std::string work_line;

  if (comment_pos != std::string::npos)
    work_line.assign(line, 0, comment_pos);
  else
    work_line = line;

  Utils::read_vector_string(work_line, number_list, ' ');

  if (number_list.size() > 0)
    ret_val = 0;  // At least one number found

  return ret_val;
  }


int load_colormap(const std::string & filename, Colormap & colormap,
		  const bool log_scale)
  {
  int ret_val = -1;  // Error state

  std::list<ValueColor> value_color_list;

  std::ifstream infile(filename.c_str());

  if (infile.good())
    {
    bool eof = false;

    std::vector<double> number_list;

    bool error = false;

    while ((!eof) && (!error))
      {
      // Read a line of data

      std::string file_line;

      std::getline(infile, file_line);

      number_list.clear();

      if (extract_numbers(file_line, number_list) >= 0)
	{
	ValueColor current_color;

	const unsigned int num_numbers = number_list.size();

	if ((num_numbers < 2) || (num_numbers > 5))
	  error = true;
	else
	  {
	  const double value = number_list[0];

	  switch (num_numbers)
	    {
	    case 2 :
	      {
	      const double gray = number_list[1];
	    
	      current_color.reset(value, gray);
	      }
	      break;

	    case 3 :
	      {
	      const double gray = number_list[1];
	      const double alpha = number_list[2];
	    
	      current_color.reset(value, gray, alpha);
	      }
	      break;

	    case 4 :
	      {
	      const double red = number_list[1];
	      const double green = number_list[2];
	      const double blue = number_list[3];

	      current_color.reset(value, red, green, blue);
	      }
	      break;

	    case 5 :
	      {
	      const double red = number_list[1];
	      const double green = number_list[2];
	      const double blue = number_list[3];

	      const double alpha = number_list[4];

	      current_color.reset(value, red, green, blue, alpha);
	      }
	      break;
	    }

	  value_color_list.push_back(current_color);
	  }
	}

      eof = infile.eof();
      }

    if (!error)
      {
      colormap.reset(value_color_list, log_scale);

      ret_val = 0;  // OK
      }
    }

  return ret_val;
  }


}  // End namespace Color
}  // End namespace QuickFlash
