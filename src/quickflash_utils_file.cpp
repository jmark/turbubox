// C++ program file quickflash_utils_file.cpp

/*
  By Nathan C. Hearn
     March 23, 2008

  Utilities for file handling.
*/


#include "quickflash_utils_file.hpp"
#include <string>


namespace QuickFlash
{
namespace Utils
{

const char directory_delimiter_unix = '/';
const char directory_delimiter_dos = '\\';


void split_name_path_file(const std::string & path_filename, 
			  std::string & path,
			  std::string & filename, 
			  const char directory_delimiter)
  {
  const std::string::size_type last_dir_marker 
    = path_filename.rfind(directory_delimiter);

  if (last_dir_marker == std::string::npos)
    {
    path = "";  // Fully-unqualified filename (no path)
    filename = path_filename;
    }
  else
    {
    const std::string::size_type str_len = path_filename.size();
    const std::string::size_type start_pos = last_dir_marker + 1;

    path = path_filename.substr(0, start_pos);

    if (start_pos < str_len)
      filename = path_filename.substr(start_pos);
    else
      filename = "";
    }
  }


}  // End namespace Utils
}  // End namespace QuickFlash
