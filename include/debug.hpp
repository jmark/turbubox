// C++ header file debug.hpp

/*
  By Nathan C. Hearn
     October 30, 2006

  Debugging-related classes.
*/


#ifndef DEBUG_HPP
#define DEBUG_HPP


#include <iostream>
using std::cerr;
using std::endl;

#include <iomanip>
using std::setprecision;


template<typename T>
std::ostream & operator<<(std::ostream & strm, const std::vector<T> & vect)
  {
  const unsigned int len = vect.size();
  
  if (len > 0)
    {
    strm << vect[0];

    for (unsigned int i = 1; i < len; i++)
      strm << " " << vect[i];
    }

  return strm;
  }



#endif  // DEBUG_HPP
