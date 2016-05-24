// C++ header file quickflash_file_datafile.hpp

/*
  By Nathan C. Hearn
     October 4, 2008

  Flash data file reader with buffered reads.
*/


#ifndef QUICKFLASH_FILE_DATAFILE_HPP
#define QUICKFLASH_FILE_DATAFILE_HPP

#include <hdf5.h>
#include <string>
#include <vector>
#include "quickflash_file_siminfo.hpp"
#include "quickflash_file_meshinfo.hpp"
#include "quickflash_file_dataset.hpp"
#include "quickflash_registry.hpp"
#include "quickflash_file_flashdefs.hpp"
#include "quickflash_particles_data.hpp"


namespace QuickFlash
{
namespace File
{

class DataFile
  {
  public :

    DataFile(const unsigned int buffer_size=0) : 
      default_buffer_size(0), data_file_name(), file_id(-1), 
      flash_code_version(FlashUnknown), flash_file_version(0), 
      sim_info(), mesh_info(), part_data(), dataset_reg() 
      { reset(buffer_size); }

    DataFile(const std::string & filename, const unsigned int buffer_size=0,
	     const bool skip_particles=false, const bool skip_mesh=false) :
      default_buffer_size(0), data_file_name(), file_id(-1),
      flash_code_version(FlashUnknown), flash_file_version(0),
      sim_info(), mesh_info(), part_data(), dataset_reg()
      { reset(filename, buffer_size, skip_particles, skip_mesh); }

    ~DataFile() { reset(); }

    void reset(const unsigned int buffer_size=0);

    void reset(const std::string & filename, const unsigned int buffer_size=0,
	       const bool skip_particles=false, const bool skip_mesh=false);

    const std::string & get_filename() const { return data_file_name; }

    FlashVersion get_flash_code_version() const { return flash_code_version; }
    unsigned int get_flash_file_version() const { return flash_file_version; }

    bool mesh_info_read() const { return mesh_info.mesh_info_read(); }

    bool mesh_valid() const { return mesh_info.mesh_valid(); }

    const SimInfo & get_sim_info() const { return sim_info; }
    const MeshInfo & get_mesh_info() const { return mesh_info; }

    bool particles_present() const { return part_data.particles_active(); }

    const Particles::PartData & get_part_data() const { return part_data; }

    const Dataset & get_dataset(const unsigned int dataset_index) const
      { return dataset_reg[dataset_index]; }

    const Dataset & get_dataset(const std::string & dataset_name,
				const unsigned int buffer_size=0) const;

    const Dataset & get_dataset(const char dataset_name[],
				const unsigned int buffer_size=0) const
      {
      const std::string dataset_name_string = dataset_name;
      return get_dataset(dataset_name_string, buffer_size);
      }

    unsigned int load_dataset(const std::string & dataset_name, 
			      const unsigned int buffer_size=0) const;

    unsigned int load_dataset(const char dataset_name[], 
			      const unsigned int buffer_size=0) const
      {
      const std::string dataset_name_string = dataset_name;
      return load_dataset(dataset_name_string, buffer_size);
      }

    void report_dataset_stats() const
      {
      const unsigned int num_datasets = dataset_reg.get_num_entries();
  
      for (unsigned int index = 0; index < num_datasets; index++)
	dataset_reg[index].report_all_stats();
      }

    void report_dataset_stats(const std::string & dataset_name) const
      { dataset_reg[dataset_name].report_all_stats(); }

  private :

    unsigned int default_buffer_size;

    std::string data_file_name;

    hid_t file_id;

    FlashVersion flash_code_version;
    unsigned int flash_file_version;

    SimInfo sim_info;
    MeshInfo mesh_info;

    Particles::PartData part_data;

    mutable Registry::Registry<std::string, Dataset> dataset_reg;
  };


}  // End namespace File
}  // End namespace QuickFlash


#endif  // QUICKFLASH_FILE_DATAFILE_HPP
