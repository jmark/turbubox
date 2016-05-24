// C++ header file quickflash_utils_option.hpp

/*
  By Nathan C. Hearn
     October 8, 2008

  Keyword-value pair container.
*/


#ifndef QUICKFLASH_UTILS_OPTION_HPP
#define QUICKFLASH_UTILS_OPTION_HPP


#include <string>
#include <sstream>
#include <vector>
#include "quickflash_utils_text.hpp"


namespace QuickFlash
{
namespace Utils
{

// Class Option

class Option
  {
  public :

    Option() : key(), value() { }

    Option(const char * option_key, const char * option_value) :
      key(), value()
      { reset(option_key, option_value); }

    Option(const std::string & option_key, 
	   const std::string & option_value) :
      key(), value()
      { reset(option_key, option_value); }

    Option(const Option & source) :
      key(), value()
      { reset(source); }

    ~Option() { }

    Option & operator=(const Option & source)
      {
      reset(source);
      return *this;
      }

    void reset()
      {
      key.clear();
      value.clear();
      }

    void reset(const char * option_key)
      { reset(option_key, 0); }

    void reset(const char * option_key, const char * option_value)
      {
      reset();

      if (option_key != 0)
	key = option_key;

      if (option_value != 0)
	value = option_value;
      }

    void reset(const std::string & option_key)
      {
      key = option_key;
      value = "";  // No value
      }

    void reset(const std::string & option_key, 
	       const std::string & option_value)
      {
      key = option_key;
      value = option_value;
      }

    void reset(const char * option_key, const std::string & option_value)
      {
      reset();

      if (option_key != 0)
	key = option_key;

      value = option_value;
      }

    void reset(const std::string & option_key, const char * option_value)
      {
      reset();

      key = option_key;

      if (option_value != 0)
	value = option_value;
      }

    void reset(const Option & source)
      {
      if (&source != this)
	{
	key = source.key;
	value = source.value;
	}
      }

    bool key_present() const { return (key != ""); }
    bool value_present() const { return (value != ""); }

    const std::string & get_key() const { return key; }
    const std::string & get_value() const { return value; }

    void get_value(std::string & option_value) const
      { option_value = value; }

    void get_value(unsigned int & option_value) const
      {
      std::stringstream sstream;

      sstream.str(value);

      sstream >> option_value;
      }

    void get_value(int & option_value) const
      {
      std::stringstream sstream;

      sstream.str(value);

      sstream >> option_value;
      }

    void get_value(double & option_value) const
      {
      std::stringstream sstream;

      sstream.str(value);

      sstream >> option_value;
      }

    void get_value(std::vector<unsigned int> & option_value) const
      { read_vector_string(value, option_value); }

    void get_value(std::vector<int> & option_value) const
      { read_vector_string(value, option_value); }

    void get_value(std::vector<double> & option_value) const
      { read_vector_string(value, option_value); }

    void get_value(std::vector<std::string> & option_value) const
      { read_vector_string(value, option_value); }

  private :

    std::string key;
    std::string value;
  };


}  // End namespace Utils
}  // End namespace QuickFlash


#endif  // QUICKFLASH_UTILS_OPTION_HPP
