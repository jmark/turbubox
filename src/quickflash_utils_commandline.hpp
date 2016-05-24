// C++ header file quickflash_utils_commandline.hpp

/*
  By Nathan C. Hearn
     February 17, 2008

  Command-line-related utilities.
*/


#ifndef QUICKFLASH_UTILS_COMMANDLINE_HPP
#define QUICKFLASH_UTILS_COMMANDLINE_HPP


#include <string>
#include <vector>
#include <list>
#include "quickflash_utils_option.hpp"
#include "quickflash_utils_text.hpp"


namespace QuickFlash
{
namespace Utils
{

extern const char command_line_option_prefix;
extern const std::string::size_type command_line_option_prefix_count;
extern const char command_line_option_delimiter;


bool is_option(const std::string & argument);

int extract_option(const std::string & argument, std::string & option_key,
		   std::string & option_value);

int extract_option(const std::string & argument, Option & option_object);

void parse_command_line(const int argc, char * const argv[], 
			std::vector<std::string> & argument_list,
			std::list<Option> & option_list,
			const bool skip_first_argv_item=true);

void parse_command_line(const int argc, char * const argv[], 
			std::vector<Option> & argument_list,
			std::list<Option> & option_list,
			const bool skip_first_argv_item=true);

    
}  // End namespace Utils
}  // End namespace QuickFlash


#endif  // QUICKFLASH_UTILS_COMMANDLINE_HPP
