// C++ header file quickflash_particles_data.hpp

/*
  By Nathan C. Hearn
     June 27, 2007

  Simulation information regarding particles.
*/


#ifndef QUICKFLASH_PARTICLES_DATA_HPP
#define QUICKFLASH_PARTICLES_DATA_HPP


#include <hdf5.h>
#include <string>
#include <vector>
#include <map>
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Particles
{

class PartData
  {
  public :

    PartData();

    PartData(const hid_t file_id);

    PartData(const PartData & source);

    ~PartData() { }

    PartData & operator=(const PartData & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const hid_t file_id);

    void reset(const PartData & source);

    bool particles_active() const { return parts_active; }

    unsigned int get_num_particles() const { return num_particles; }

    bool particle_present(const unsigned int particle_tag) const;

    unsigned int get_particle_index(const unsigned int particle_tag) const
      {
      const std::map<unsigned int, unsigned int>::const_iterator tag_iter
	= particle_tag_index_map.find(particle_tag);

      if (tag_iter == particle_tag_index_map.end())
	throw Except("Particle variable not present", __FILE__, __LINE__);

      return tag_iter->second;
      }
      
    unsigned int get_particle_tag(const unsigned int particle_index) const
      { return particle_tags[particle_index]; }

    const std::vector<unsigned int> & get_particle_tag_list() const
      { return particle_tags; }

    unsigned int get_num_variables() const { return num_variables; }

    bool variable_present(const std::string & var_name) const;

    unsigned int get_variable_index(const std::string & var_name) const;

    const std::string & get_variable_name(const unsigned int var_index) const
      { return variable_names[var_index]; }

    const std::vector<std::string> & get_variable_names() const
      { return variable_names; }

    void get_variable_data(const unsigned int var_index, 
			   std::vector<double> & data) const;

    void get_variable_data(const std::string & var_name, 
			   std::vector<double> & data) const;

    void get_particle_data_index(const unsigned int part_index,
				 std::vector<double> & data) const;

    void get_particle_data_tag(const unsigned int part_tag,
			       std::vector<double> & data) const
      {
      const unsigned int part_index = get_particle_index(part_tag);

      get_particle_data_index(part_index, data);
      }


  private :

    void read_particle_tags(const hid_t file_id);

  private :

    bool parts_active;

    hid_t part_file_id;  // Just a reference -- DO NOT CLOSE!!!

    unsigned int num_particles;

    std::vector<unsigned int> particle_tags;
    std::map<unsigned int, unsigned int> particle_tag_index_map;

    unsigned int num_variables;

    std::vector<std::string> variable_names;
    std::map<std::string, unsigned int> variable_name_index_map;
  };

}  // End namespace Particles
}  // End namespace QuickFlash


#endif  // QUICKFLASH_PARTICLES_DATA_HPP
