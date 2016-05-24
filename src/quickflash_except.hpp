// C++ header file quickflash_except.hpp

/*
  By Nathan C. Hearn
     August 23, 2006

  Simple exception class.
*/


#ifndef QUICKFLASH_EXCEPT_HPP
#define QUICKFLASH_EXCEPT_HPP


#include <exception>
#include <string>
#include <iostream>
#include <sstream>
#include <list>
#include <utility>


// #define EXCEPT(MESSAGE_STR) (Except(MESSAGE_STR, __FILE__, __LINE__))


namespace QuickFlash
{

class Except : public std::exception
  {
  public :

    Except(const std::string & message, const char source_file[]="", 
	   const int line_number=-1, const bool skip_announce=false) :
      exception(), already_announced(false), msg(message), 
      filename(source_file), line(line_number)
      { 
      if (!skip_announce)
	announce(); 
      }

    Except(const char message[], const char source_file[]="", 
	   const int line_number=-1, const bool skip_announce=false) :
      exception(), already_announced(false), msg(message), 
      filename(source_file), line(line_number)
      {
      if (!skip_announce)
	announce(); 
      }

    Except(const Except & source) :
      exception(), already_announced(source.already_announced), 
      msg(source.msg), filename(source.filename), line(source.line)
      {
      if (!already_announced)
	announce();
      }

    virtual ~Except() throw () { }

    Except & operator=(const Except & source)
      {
      if (&source != this)
	{
	already_announced = source.already_announced;

	msg = source.msg;

	filename = source.filename;
	line = source.line;
	}

      return *this;
      }

    void add_keyword_value(const std::string & keyword)
      {
      const std::pair<std::string, std::string> keyword_value(keyword, "");

      keyword_value_list.push_back(keyword_value);
      }

    template<typename T>
    void add_keyword_value(const std::string & keyword, const T & value)
      {
      std::stringstream value_stream;

      value_stream << value;

      const std::pair<std::string, std::string> 
	keyword_value(keyword, value_stream.str());

      keyword_value_list.push_back(keyword_value);
      }

    const std::string & get_message() const { return msg; }

    const std::string & get_source_filename() const { return filename; }

    int get_line_number() const { return line; }

    void announce() const
      {
      std::stringstream sstream;

      sstream << "Error [ " << msg << " ]" << std::endl;

      if (filename.size() > 0)
	{
	sstream << "  - Source file [ " << filename << " ]";

	if (line >= 0)
	  sstream << " line [ " << line << " ]";
	}

      std::cerr << sstream.str() << std::endl;

      std::list< std::pair<std::string, std::string> >::const_iterator 
	keyword_value_iter = keyword_value_list.begin();

      while (keyword_value_iter != keyword_value_list.end())
	{
	const std::string & keyword = keyword_value_iter->first;
	const std::string & value = keyword_value_iter->second;

	std::cerr << "  - " << keyword << " " << value << std::endl;

	++keyword_value_iter;
	}

      std::cerr << std::endl;

      already_announced = true;
      }

  private :

    mutable bool already_announced;

    std::string msg;

    std::string filename;
    int line;

    std::list< std::pair<std::string, std::string > > keyword_value_list;
  };


}  // End namespace QuickFlash


#endif  // QUICKFLASH_EXCEPT_HPP
