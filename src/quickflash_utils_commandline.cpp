// C++ program file quickflash_utils_commandline.hpp

/*
  By Nathan C. Hearn
     February 17, 2008

  Command-line-related utilities.
*/


#include "quickflash_utils_commandline.hpp"
#include <string>
#include <vector>
#include <list>
#include "quickflash_utils_option.hpp"
#include "quickflash_utils_text.hpp"
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Utils
{

const char command_line_option_prefix = '-';
const std::string::size_type command_line_option_prefix_count = 2;
const char command_line_option_delimiter = '=';


bool is_option(const std::string & argument)
  {
  bool ret_val = false;

  // For right now, assume no characters before prefix

  if (argument.size() >= command_line_option_prefix_count)
    {
    ret_val = true;  // Could be an option ...

    for (std::string::size_type i = 0; i < command_line_option_prefix_count;
	 i++)
      if (argument[i] != command_line_option_prefix)
	ret_val = false;
    }

  return ret_val;
  }
  

int extract_option(const std::string & argument, std::string & option_key,
		   std::string & option_value)
  {
  int ret_val = -1;  // Error state

  if (is_option(argument))
    {
    const std::string::size_type arg_len = argument.size();

    if (arg_len > command_line_option_prefix_count)
      {
      // See if there is a delimiter

      const std::string::size_type delimiter_pos 
	= argument.find(command_line_option_delimiter, 0);

      // Find the end of the key and copy
      const std::string::size_type key_end 
	= (delimiter_pos == std::string::npos) ? arg_len : delimiter_pos;

      const std::string::size_type key_len 
	= key_end - command_line_option_prefix_count;

      option_key = argument.substr(command_line_option_prefix_count, key_len);

      // Find the value, if there is one

      if (key_end < (arg_len - 1))
	{
	const std::string::size_type value_begin = key_end + 1;
	const std::string::size_type value_len = arg_len - value_begin;

	option_value = argument.substr(value_begin, value_len);
	}
      else
	option_value = "";
      }
    else
      {
      option_key = "";
      option_value = "";
      }

    ret_val = 0;  // OK
    }

  return ret_val;
  }


int extract_option(const std::string & argument, Option & option_object)
  {
  std::string key;
  std::string value;

  const int ret_val = extract_option(argument, key, value);

  if (ret_val >= 0)
    option_object.reset(key, value);

  return ret_val;
  }


void parse_command_line(const int argc, char * const argv[], 
			std::vector<std::string> & argument_list,
			std::list<Option> & option_list,
			const bool skip_first_argv_item)
  {
  // Parse argv, collecting all non-options into argument_list, and all
  // options into option_list
  //
  // Note: All "options" after an empty-key option (i.e., nothing following
  // "--") are considered arguments

  if (argv == 0)
    throw Except("Argument list missing", __FILE__, __LINE__);

  argument_list.clear();
  option_list.clear();

  int argptr = (skip_first_argv_item) ? 1 : 0;

  bool only_args_remain = false;

  while (argptr < argc)
    {
    const char * arg_item = argv[argptr++];

    if (arg_item == 0)
      throw Except("Argument item missing", __FILE__, __LINE__);

    // OTHER ERROR CHECKING???

    const std::string arg_string = arg_item;

    if ((only_args_remain) || (!is_option(arg_string)))
      argument_list.push_back(arg_string);
    else
      {
      Option new_option;

      if (extract_option(arg_string, new_option) < 0)
	throw Except("Unable to parse option", __FILE__, __LINE__);

      if (!new_option.key_present())
	only_args_remain = true;
      else
	option_list.push_back(new_option);
      }
    }
  }


void parse_command_line(const int argc, char * const argv[], 
			std::vector<Option> & argument_list,
			std::list<Option> & option_list,
			const bool skip_first_argv_item)
  {
  std::vector<std::string> argument_string_list;

  parse_command_line(argc, argv, argument_string_list, option_list,
		     skip_first_argv_item);

  const unsigned int num_args = argument_string_list.size();

  argument_list.resize(num_args);

  for (unsigned int index = 0; index < num_args; index++)
    argument_list[index].reset(0, argument_string_list[index]);
  }

    
}  // End namespace Utils
}  // End namespace QuickFlash
