// C++ program file quickflash_particles_data.cpp

/*
  By Nathan C. Hearn
     June 27, 2007

  Simulation information regarding particles.
*/


#include "quickflash_particles_data.hpp"
#include <hdf5.h>
#include <string>
#include <vector>
#include <map>
#include "quickflash_hdf5.hpp"
#include "quickflash_file_flashdefs.hpp"
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace Particles
{

PartData::PartData() :
  parts_active(false), part_file_id(-1), num_particles(0), particle_tags(), 
  particle_tag_index_map(), num_variables(0), variable_names(),
  variable_name_index_map()
  { }


PartData::PartData(const hid_t file_id) :
  parts_active(false), part_file_id(-1), num_particles(0), particle_tags(), 
  particle_tag_index_map(), num_variables(0), variable_names(),
  variable_name_index_map()
  { reset(file_id); }


PartData::PartData(const PartData & source) :
  parts_active(false), part_file_id(-1), num_particles(0), particle_tags(), 
  particle_tag_index_map(), num_variables(0), variable_names(),
  variable_name_index_map()
  { reset(source); }


void PartData::reset()
  {
  parts_active = false;

  part_file_id = -1;

  num_particles = 0;

  particle_tags.clear();
  particle_tag_index_map.clear();

  num_variables = 0;

  variable_names.clear();
  variable_name_index_map.clear();
  }


void PartData::reset(const hid_t file_id)
  {
  reset();

  if (file_id <= 0)
    throw Except("Improper HDF5 container ID", __FILE__, __LINE__);

  File::FlashVersion flash_version = File::Flash2;
  unsigned int file_version = 0;

  if (get_flash_version(file_id, file_version, flash_version) < 0)
    throw Except("Unable to read Flash file version", __FILE__, __LINE__);

  if (flash_version != File::Flash3)
    throw Except("Particles supported only in Flash3", __FILE__, __LINE__);

  // Search the file for the necessary components -- ASSUME FLASH3 FOR NOW

  if (HDF5::member_present(file_id, File::Particle_Variables_Name_Flash3))
    {
    if (HDF5::read_string_set(file_id, File::Particle_Variables_Name_Flash3,
			      variable_names) < 0)
      throw Except("Unable to read particle variables name list", __FILE__,
		   __LINE__);

    num_variables = variable_names.size();

    for (unsigned int index = 0; index < num_variables; index++)
      variable_name_index_map[variable_names[index]] = index;

    try
      {
      read_particle_tags(file_id);
      }
    catch (Except ex)
      {
      reset();
      throw Except(ex);
      }

    parts_active = true;

    part_file_id = file_id;  // Just a reference ... DO NOT CLOSE!!!
    }
  }


void PartData::reset(const PartData & source)
  {
  if (&source != this)
    {
    parts_active = source.parts_active;

    part_file_id = source.part_file_id;

    num_particles = source.num_particles;

    particle_tags = source.particle_tags;
    particle_tag_index_map = source.particle_tag_index_map;

    num_variables = source.num_variables;

    variable_names = source.variable_names;
    variable_name_index_map = source.variable_name_index_map;
    }
  }


bool PartData::particle_present(const unsigned int particle_tag) const
  {
  bool found_it = false;

  const std::map<unsigned int, unsigned int>::const_iterator part_iter
    = particle_tag_index_map.find(particle_tag);

  if (part_iter != particle_tag_index_map.end())
    found_it = true;

  return found_it;
  }


bool PartData::variable_present(const std::string & var_name) const
  {
  bool found_it = false;

  const std::map<std::string, unsigned int>::const_iterator var_iter
    = variable_name_index_map.find(var_name);

  if (var_iter != variable_name_index_map.end())
    found_it = true;

  return found_it;
  }


unsigned int PartData::get_variable_index(const std::string & var_name) const
  {
  const std::map<std::string, unsigned int>::const_iterator var_iter
    = variable_name_index_map.find(var_name);

  if (var_iter == variable_name_index_map.end())
    throw Except("Particle variable not present", __FILE__, __LINE__);

  return var_iter->second;
  }


void PartData::get_variable_data(const unsigned int var_index, 
				 std::vector<double> & data) const
  {
  if (!parts_active)
    throw Except("Particles not active", __FILE__, __LINE__);

  if (var_index >= num_variables)
    throw Except("Variable index out of range", __FILE__, __LINE__);

  std::vector<unsigned int> select_coords(2);

  select_coords[0] = 0;
  select_coords[1] = var_index;

  if (HDF5::read_data_stripe(part_file_id, File::Particle_Data_Name, 0,
			     select_coords, data) < 0)
    throw Except("Unable to read particle data", __FILE__, __LINE__);
  }


void PartData::get_variable_data(const std::string & var_name, 
				 std::vector<double> & data) const
  {
  const unsigned int var_index = get_variable_index(var_name);

  get_variable_data(var_index, data);
  }


void PartData::get_particle_data_index(const unsigned int part_index, 
				       std::vector<double> & data) const
  {
  if (!parts_active)
    throw Except("Particles not active", __FILE__, __LINE__);

  if (part_index >= num_particles)
    throw Except("Particle index out of range", __FILE__, __LINE__);

  std::vector<unsigned int> select_coords(2);

  select_coords[0] = part_index;
  select_coords[1] = 0;

  if (HDF5::read_data_stripe(part_file_id, File::Particle_Data_Name, 1,
			     select_coords, data) < 0)
    throw Except("Unable to read particle data", __FILE__, __LINE__);
  }


void PartData::read_particle_tags(const hid_t file_id)
  {
  // Call ONLY after variable names have been read!!

  const std::string tag_string = File::Particle_Variables_Tag_Name;

  const unsigned int particle_tag_index = get_variable_index(tag_string);

  // Build selection vector -- remember column-major storage order!!!

  std::vector<unsigned int> select_coords(2);

  select_coords[0] = 0;
  select_coords[1] = particle_tag_index;

  std::vector<double> double_tags;

  if (HDF5::read_data_stripe(file_id, File::Particle_Data_Name, 0,
			     select_coords, double_tags) < 0)
    throw Except("Unable to read particle tags", __FILE__, __LINE__);

  num_particles = double_tags.size();

  particle_tags.resize(num_particles);

  particle_tag_index_map.clear();

  for (unsigned int index = 0; index < num_particles; index++)
    {
    const unsigned int tag = static_cast<unsigned int>(double_tags[index]);

    particle_tags[index] = tag;

    if (particle_tag_index_map.find(tag) != particle_tag_index_map.end())
      {
	std::cerr << "Tag " << tag << " at index " << index << " already stored at index " << particle_tag_index_map[tag] << std::endl;
      throw Except("Duplicate particle tag found", __FILE__, __LINE__);
      }

    particle_tag_index_map[tag] = index;
    }
  }


}  // End namespace Particles
}  // End namespace QuickFlash
