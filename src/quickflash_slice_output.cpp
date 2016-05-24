// C++ program file quickflash_slice_output.cpp

/*
  By Nathan C. Hearn
     April 23, 2008

  File output routines (ASCII) for SliceGrid data.
*/


#include "quickflash_slice_output.hpp"
#include <string>
#include <list>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "quickflash_slice_slicegrid.hpp"
#include "quickflash_utils_text.hpp"


namespace QuickFlash
{
namespace Slice
{

const double rad_to_deg = 180.0 / (M_PI);


void write_grid(const std::string & filename, const SliceGrid & data, 
		const bool add_comment_marker, const char comment_marker)
  {
  std::list<std::string> empty_metadata_list;

  write_grid(filename, data, empty_metadata_list, add_comment_marker, 
	     comment_marker);
  }


void write_grid(const std::string & filename, 
		const SliceGrid & data,
		const std::list<std::string> & metadata_list,
		const bool add_comment_marker, const char comment_marker)
  {
  using std::endl;

  std::string comment_prefix = "";

  if (add_comment_marker)
    comment_prefix += comment_marker + " ";

  std::ofstream outfile(filename.c_str());

  outfile << comment_prefix << "GridDataFile " << filename << endl;

  std::list<std::string>::const_iterator str_iter = metadata_list.begin();
  std::list<std::string>::const_iterator end_str_iter = metadata_list.end();

  while (str_iter != end_str_iter)
    {
    outfile << comment_prefix << *str_iter << endl;
    ++str_iter;
    }

  const std::vector<double> & grid_center = data.get_grid_center();

  outfile << comment_prefix << "GridCenter " << grid_center << endl;

  const std::vector<unsigned int> & grid_dims = data.get_grid_dims();

  const unsigned int num_grid_dims = grid_dims.size();

  if (num_grid_dims != 2)
    throw Except("Incompatible grid dimensions", __FILE__, __LINE__);

  outfile << comment_prefix << "GridPixels " << grid_dims << endl;

  const std::vector<double> & grid_width = data.get_slice_grid_width();

  outfile << comment_prefix << "GridWidth " << grid_width << endl;

  const double roll_angle_deg = data.get_roll_angle() * rad_to_deg;
  const double polar_angle_deg = data.get_polar_angle() * rad_to_deg;
  const double azim_angle_deg = data.get_azimuthal_angle() * rad_to_deg;

  outfile << comment_prefix << "RollAngle " << roll_angle_deg << endl;
  outfile << comment_prefix << "PolarAngle " << polar_angle_deg << endl;
  outfile << comment_prefix << "AzimuthalAngle " << azim_angle_deg << endl;

  double min_value = 0.0;
  double max_value = 0.0;

  data.get_minmax_values(min_value, max_value);

  outfile << comment_prefix << "MinValue " << min_value << endl;
  outfile << comment_prefix << "MaxValue " << max_value << endl;

  // Write out the data

  outfile << comment_prefix << "SliceData" << endl;

  const unsigned int max_i = grid_dims[0];
  const unsigned int max_j = grid_dims[1];

  std::vector<unsigned int> pixel_index(2);

  for (unsigned int i = 0; i < max_i; i++)
    {
    pixel_index[0] = i;

    for (unsigned int j = 0; j < max_j; j++)
      {
      pixel_index[1] = j;

      outfile << data[pixel_index] << " ";
      }

    outfile << endl;
    }

  outfile << endl;
  }


}  // End namespace Slice
}  // End namespace QuickFlash
