// C++ header file quickflash_utils_text.hpp

/*
  By Nathan C. Hearn
     February 11, 2008

  Text-related utilities for QuickFlash.
*/


#ifndef QUICKFLASH_UTILS_TEXT_HPP
#define QUICKFLASH_UTILS_TEXT_HPP


#include <iostream>
#include <vector>


namespace QuickFlash
{
namespace Utils
{

extern const char vector_element_separator;


// Vector-related functions 

template<class T>
void get_vector_string(const std::vector<T> & vect, std::string & vect_string,
		       const char element_separator=vector_element_separator);

void get_vector_string(const std::vector<bool> & vect, 
		       std::string & vect_string,
		       const char element_separator=vector_element_separator);

template<class T>
void read_vector_string(const std::string & vect_string, 
			std::vector<T> & vect,
                        const char element_separator=vector_element_separator);

void read_vector_string(const std::string & vect_string, 
			std::vector<bool> & vect,
                        const char element_separator=vector_element_separator);


}  // End namespace Utils
}  // End namespace QuickFlash


namespace std
{

// Operators for use with iostream objects

template<class T>
ostream & operator<<(ostream & ostr, const vector<T> & vect);

template<class T>
istream & operator>>(istream & istr, vector<T> & vect);


}  // End namespace std


// Include the function definitions

#include "quickflash_utils_text.tcc"


#endif  // QUICKFLASH_UTILS_TEXT_HPP

