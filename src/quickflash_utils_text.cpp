// C++ program file quickflash_utils_text.cpp

/*
  By Nathan C. Hearn
     February 11, 2008

  Text-related utilities for QuickFlash.
*/


#include "quickflash_utils_text.hpp"
#include <vector>
#include <string>
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Utils
{

const char vector_element_separator = ',';

void get_vector_string(const std::vector<bool> & vect, 
		       std::string & vect_string, const char element_separator)
  {
  const unsigned int num_elems = vect.size();

  if (num_elems > 0)
    {
    if (vect[0])
      vect_string = "true";
    else
      vect_string = "false";

    for (unsigned int i = 1; i < num_elems; i++)
      {
      vect_string += element_separator;

      if (vect[i])
	vect_string += "true";
      else
	vect_string += "false";
      }
    }
  }


void read_vector_string(const std::string & vect_string, 
			std::vector<bool> & vect, const char element_separator)
  {
  const std::string::size_type num_chars = vect_string.size();

  // string_end points to the location AFTER the last character in vect_string

  std::string::size_type string_begin = 0;
  std::string::size_type string_end = num_chars;

  // Ignore any space characters at the beginning

  bool found_begin = false;

  while ((!found_begin) && (string_begin < string_end))
    {
    if (vect_string[string_begin] != ' ')
      found_begin = true;
    else
      string_begin++;
    }

  // Ignore any space characters at the end

  bool found_end = false;

  while (!found_end)
    {
    if (string_end == string_begin)
      found_end = true;
    else
      {
      const std::string::size_type prev_char_index = string_end - 1;

      if (vect_string[prev_char_index] != ' ')
	found_end = true;
      else
	string_end = prev_char_index;
      }
    }

  // Read the vector elements

  vect.clear();

  if (string_end > string_begin)
    {
    std::string::size_type current_string_index = string_begin;

    bool finished = false;

    while (!finished)
      {
      // Find the next separator

      const std::string::size_type next_separator 
	= vect_string.find(element_separator, current_string_index);

      // current_end_index points either to the next separator or string_end

      std::string::size_type current_end_index
	= (next_separator == std::string::npos) ? string_end : next_separator;

      const std::string::size_type current_length 
	= current_end_index - current_string_index;

      const std::string elem_str = vect_string.substr(current_string_index, 
						      current_length);

      bool data_elem = false;

      if ((elem_str == "true") || (elem_str == "True") || (elem_str == "T")
	  || (elem_str == "t"))
	data_elem = true;
      else if ((elem_str == "false") || (elem_str == "False") 
	       || (elem_str == "F") || (elem_str == "f"))
	data_elem = false;
      else
	Except("Unable to interpret boolean text string", __FILE__, __LINE__);

      vect.push_back(data_elem);

      // Move to the last consecutive space

      bool found_last_space = false;

      while (!found_last_space)
	{
	const std::string::size_type next_index = current_end_index + 1;

	if (next_index < string_end)
	  {
	  const char next_character = vect_string[next_index];

	  if (next_character != ' ')
	    found_last_space = true;
	  else
	    current_end_index = next_index;
	  }
	else
	  found_last_space = true;
	}

      current_string_index = current_end_index + 1;  // Skip the separator

      if (current_string_index >= string_end)
	finished = true;
      }
    }
  }


}  // End namespace Utils
}  // End namespace QuickFlash
