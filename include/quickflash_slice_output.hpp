// C++ header file quickflash_slice_output.hpp

/*
  By Nathan C. Hearn
     April 23, 2008

  File output routines (ASCII) for SliceGrid data.
*/


#ifndef QUICKFLASH_SLICE_OUTPUT_HPP
#define QUICKFLASH_SLICE_OUTPUT_HPP


#include <string>
#include <list>
#include "quickflash_slice_slicegrid.hpp"


namespace QuickFlash
{
namespace Slice
{

extern const double rad_to_deg;


void write_grid(const std::string & filename, const SliceGrid & data, 
		const bool add_comment_marker=false, 
		const char comment_marker='#');

void write_grid(const std::string & filename, 
		const SliceGrid & data,
		const std::list<std::string> & metadata_list,
		const bool add_comment_marker=false,
		const char comment_marker='#');


}  // End namespace Slice
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SLICE_OUTPUT_HPP
