// C++ program file quickflash_datafile.cpp

/*
  By Nathan C. Hearn
     September 2, 2006

  Flash data file reader with buffered reads.
*/


#include "quickflash_file_datafile.hpp"
#include <hdf5.h>
#include <string>
#include "quickflash_file_flashdefs.hpp"
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace File
{

void DataFile::reset(const unsigned int buffer_size)
  {
  if (buffer_size > 0)
    default_buffer_size = buffer_size;
  else if (default_buffer_size < 1)
    default_buffer_size = 1;

  data_file_name = "";

  dataset_reg.reset();  // Reset the registry before closing the file

  if (file_id >= 0)
    {
    H5Fclose(file_id);
    file_id = -1;
    }

  flash_code_version = FlashUnknown;
  flash_file_version = 0;

  sim_info.reset();
  mesh_info.reset();

  part_data.reset();
  }


void DataFile::reset(const std::string & filename,
		     const unsigned int buffer_size,
		     const bool skip_particles, const bool skip_mesh)
  {
  reset();

  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0)
    {
    Except exception("Unable to open file", __FILE__, __LINE__, true);

    exception.add_keyword_value("Filename", filename);
    exception.announce();

    throw exception;
    }

  if (get_flash_version(file_id, flash_file_version, flash_code_version) < 0)
    {
    reset();
    throw Except("Unable to get Flash file version", __FILE__, __LINE__);
    }

  sim_info.reset(file_id, flash_code_version, flash_file_version, skip_mesh);

  if (!skip_mesh)
    {
    try
      {
      mesh_info.reset(file_id, sim_info);
      }
    catch (Except ex)
      {
      // Just run without mesh info for now ...

      mesh_info.reset();
      }
    }

  if (!skip_particles)
    {
    try
      {
      part_data.reset(file_id);
      }
    catch (Except ex)
      {
      // Just run without particles for now

      part_data.reset();
      }
    }

  if (buffer_size > 0)
    default_buffer_size = buffer_size;
  else if (default_buffer_size < 1)
    default_buffer_size = 1;

  data_file_name = filename;
  }


const Dataset & DataFile::get_dataset(const std::string & dataset_name,
				      const unsigned int buffer_size) const
  {
  if (!mesh_valid())
    throw Except("Mesh information not available", __FILE__, __LINE__);

  std::string padded_dataset_name;

  pad_variable_name_file_version(flash_file_version, dataset_name, 
				 padded_dataset_name);

  if (!(dataset_reg.key_present(padded_dataset_name)))
    {
    unsigned int current_buffer_size = buffer_size;

    if (current_buffer_size < 1)
      current_buffer_size = default_buffer_size;

    load_dataset(padded_dataset_name, current_buffer_size);
    }

  return dataset_reg[padded_dataset_name];
  }


unsigned int DataFile::load_dataset(const std::string & dataset_name, 
				    const unsigned int buffer_size) const
  {
  if (!mesh_valid())
    throw Except("Mesh information not available", __FILE__, __LINE__);

  const unsigned int space_dims = mesh_info.get_dims();

  std::string padded_dataset_name;

  pad_variable_name_file_version(flash_file_version, dataset_name, 
				 padded_dataset_name);

  unsigned int current_buffer_size = buffer_size;

  if (current_buffer_size < 1)
    current_buffer_size = default_buffer_size;

  const unsigned int dataset_index = dataset_reg.add_new(padded_dataset_name);

  Dataset & dset = dataset_reg[dataset_index];

  dset.reset(file_id);
  dset.open_dataset(padded_dataset_name, space_dims, current_buffer_size);

  return dataset_index;
  }


}  // End namespace File
}  // End namespace QuickFlash
