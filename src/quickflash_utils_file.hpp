// C++ header file quickflash_utils_file.hpp

/*
  By Nathan C. Hearn
     March 23, 2008

  Utilities for file handling.
*/


#ifndef QUICKFLASH_UTILS_FILE_HPP
#define QUICKFLASH_UTILS_FILE_HPP


#include <string>


namespace QuickFlash
{
namespace Utils
{

extern const char directory_delimiter_unix;
extern const char directory_delimiter_dos;

void split_name_path_file(const std::string & path_filename, 
			  std::string & path,
			  std::string & filename, 
			  const char 
			    directory_delimiter=directory_delimiter_unix);


}  // End namespace Utils
}  // End namespace QuickFlash


#endif  // QUICKFLASH_UTILS_FILE_HPP
