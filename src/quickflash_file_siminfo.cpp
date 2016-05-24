// C++ program file quickflash_file_siminfo.cpp

/*
  By Nathan C. Hearn
     October 4, 2008

  Container for general simulation information.
*/


#include "quickflash_file_siminfo.hpp"
#include <hdf5.h>
#include <string>
#include <map>
#include "quickflash_file_flashdefs.hpp"
#include "quickflash_geometry.hpp"
#include "quickflash_except.hpp"
#include "quickflash_hdf5.hpp"


namespace QuickFlash
{
namespace File
{

// Class SimInfo

SimInfo::SimInfo() :
  flash_version(FlashUnknown), file_version(0), valid_mesh(false), 
  valid_particles(false), mesh_dims(0), total_blocks(0), total_particles(0),
  sim_time(0.0), sim_dt(0.0), step_count(0), block_dims(), block_cell_count(0),
  num_mesh_variables(0), mesh_variable_names(), base_block_dims(), 
  volume_min_bounds(), volume_max_bounds(), volume_center(), volume_width(),
  logical_runtime_parameters(), integer_runtime_parameters(),
  real_runtime_parameters(), string_runtime_parameters(),
  geom_type(Geometry::Cartesian)
  { }


SimInfo::SimInfo(const hid_t file_id, const FlashVersion flash_code_version,
		 const unsigned int flash_file_version, const bool skip_mesh) :
  flash_version(FlashUnknown), file_version(0), valid_mesh(false), 
  valid_particles(false), mesh_dims(0), total_blocks(0), total_particles(0),
  sim_time(0.0), sim_dt(0.0), step_count(0), block_dims(), block_cell_count(0),
  num_mesh_variables(0), mesh_variable_names(), base_block_dims(), 
  volume_min_bounds(), volume_max_bounds(), volume_center(), volume_width(),
  logical_runtime_parameters(), integer_runtime_parameters(),
  real_runtime_parameters(), string_runtime_parameters(),
  geom_type(Geometry::Cartesian)
  { reset(file_id, flash_code_version, flash_file_version, skip_mesh); }


SimInfo::SimInfo(const SimInfo & source) :
  flash_version(FlashUnknown), file_version(0), valid_mesh(false), 
  valid_particles(false), mesh_dims(0), total_blocks(0), total_particles(0),
  sim_time(0.0), sim_dt(0.0), step_count(0), block_dims(), block_cell_count(0),
  num_mesh_variables(0), mesh_variable_names(), base_block_dims(),
  volume_min_bounds(), volume_max_bounds(), volume_center(), volume_width(),
  logical_runtime_parameters(), integer_runtime_parameters(),
  real_runtime_parameters(), string_runtime_parameters(),
  geom_type(Geometry::Cartesian)
  { reset(source); }


void SimInfo::reset()
  {
  flash_version = FlashUnknown;
  file_version = 0;

  valid_mesh = false;
  valid_particles = false;

  mesh_dims = 0;

  total_blocks = 0;
  total_particles = 0;

  sim_time = 0.0;
  sim_dt = 0.0;

  step_count = 0;

  block_dims.clear();
  block_cell_count = 0;

  num_mesh_variables = 0;
  mesh_variable_names.clear();

  base_block_dims.clear();

  volume_min_bounds.clear();
  volume_max_bounds.clear();

  volume_center.clear();
  volume_width.clear();

  logical_runtime_parameters.clear();
  integer_runtime_parameters.clear();
  real_runtime_parameters.clear();
  string_runtime_parameters.clear();

  geom_type = Geometry::Cartesian;
  }


void SimInfo::reset(const hid_t file_id, const FlashVersion flash_code_version,
		    const unsigned int flash_file_version, 
		    const bool skip_mesh)
  {
  reset();

  if (file_id < 0)
    throw Except("Invalid HDF5 file handle", __FILE__, __LINE__);

  flash_version = flash_code_version;
  file_version = flash_file_version;

  if (skip_mesh)
    {
    valid_mesh = false;
    mesh_dims = 0;
    }
  else if (get_mesh_dimensionality(file_id, file_version, mesh_dims) < 0)
    {
    valid_mesh = false;
    mesh_dims = 0;
    }
  else
    valid_mesh = true;

  read_simulation_params(file_id);

  if (total_particles > 0)
    valid_particles = true;
  else
    valid_particles = false;

  read_runtime_params(file_id);

  // Read the mesh variable names

  if (valid_mesh)
    {
    if (HDF5::read_string_set(file_id, Variables_Name, 
			      mesh_variable_names) < 0)
      throw Except("Unable to read mesh variable names", __FILE__, __LINE__);

    num_mesh_variables = mesh_variable_names.size();
    }
  else
    {
    mesh_variable_names.clear();
    num_mesh_variables = 0;
    }

  // Determine the physical geometry

  // const std::string geometry_string = Geometry_Label;

  if (!(string_runtime_param_present(Geometry_Label)))
    throw Except("Geometry description string not present", __FILE__, 
		 __LINE__);

  const std::string & geometry_type 
    = get_string_runtime_param(Geometry_Label);

  if (Geometry::get_geometry_type(geometry_type, geom_type) < 0)
    throw Except("Geometry type not recognized", __FILE__, __LINE__);
  }


void SimInfo::reset(const SimInfo & source)
  {
  if (&source != this)
    {
    flash_version = source.flash_version;
    file_version = source.file_version;

    valid_mesh = source.valid_mesh;
    valid_particles = source.valid_particles;

    mesh_dims = source.mesh_dims;

    total_blocks = source.total_blocks;
    total_particles = source.total_particles;

    sim_time = source.sim_time;
    sim_dt = source.sim_dt;

    step_count = source.step_count;

    block_dims = source.block_dims;
    block_cell_count = source.block_cell_count;

    num_mesh_variables = source.num_mesh_variables;
    mesh_variable_names = source.mesh_variable_names;

    base_block_dims = source.base_block_dims;

    volume_min_bounds = source.volume_min_bounds;
    volume_max_bounds = source.volume_max_bounds;

    volume_center = source.volume_center;
    volume_width = source.volume_width;

    logical_runtime_parameters = source.logical_runtime_parameters;
    integer_runtime_parameters = source.integer_runtime_parameters;
    real_runtime_parameters = source.real_runtime_parameters;
    string_runtime_parameters = source.string_runtime_parameters;

    geom_type = source.geom_type;
    }
  }


bool SimInfo::logical_runtime_param_present(const std::string & name) const
  {
  bool ret_val = false;

  const std::map<std::string, bool>::const_iterator iter 
      = logical_runtime_parameters.find(name);

  if (iter != logical_runtime_parameters.end())
    ret_val = true;

  return ret_val;
  }


bool SimInfo::get_logical_runtime_param(const std::string & name) const
  {
  const std::map<std::string, bool>::const_iterator iter
    = logical_runtime_parameters.find(name);

  if (iter == logical_runtime_parameters.end())
    throw Except("Boolean runtime parameter not present", __FILE__, __LINE__);

  return iter->second;
  }
  

bool SimInfo::integer_runtime_param_present(const std::string & name) const
  {
  bool ret_val = false;

  const std::map<std::string, int>::const_iterator iter 
      = integer_runtime_parameters.find(name);

  if (iter != integer_runtime_parameters.end())
    ret_val = true;

  return ret_val;
  }


int SimInfo::get_integer_runtime_param(const std::string & name) const
  {
  const std::map<std::string, int>::const_iterator iter
    = integer_runtime_parameters.find(name);

  if (iter == integer_runtime_parameters.end())
    throw Except("Integer runtime parameter not present", __FILE__, __LINE__);

  return iter->second;
  }
  

bool SimInfo::real_runtime_param_present(const std::string & name) const
  {
  bool ret_val = false;

  const std::map<std::string, double>::const_iterator iter 
      = real_runtime_parameters.find(name);

  if (iter != real_runtime_parameters.end())
    ret_val = true;

  return ret_val;
  }


double SimInfo::get_real_runtime_param(const std::string & name) const
  {
  const std::map<std::string, double>::const_iterator iter
    = real_runtime_parameters.find(name);

  if (iter == real_runtime_parameters.end())
    throw Except("Floating point runtime parameter not present", __FILE__, 
		 __LINE__);

  return iter->second;
  }


bool SimInfo::string_runtime_param_present(const std::string & name) const
  {
  bool ret_val = false;

  const std::map<std::string, std::string>::const_iterator iter 
      = string_runtime_parameters.find(name);

  if (iter != string_runtime_parameters.end())
    ret_val = true;

  return ret_val;
  }


const std::string & SimInfo::get_string_runtime_param(const std::string 
							 & name) const
  {
  const std::map<std::string, std::string>::const_iterator iter
    = string_runtime_parameters.find(name);

  if (iter == string_runtime_parameters.end())
    throw Except("Floating point runtime parameter not present", __FILE__, 
		 __LINE__);

  return iter->second;
  }


bool SimInfo::mesh_variable_present(const char variable_name[]) const
  {
  bool ret_val = false;

  const unsigned int num_vars = mesh_variable_names.size();

  unsigned int index = 0;

  while ((!ret_val) && (index < num_vars))
    {
    if (mesh_variable_names[index] == variable_name)
      ret_val = true;
    else
      index++;
    }

  return ret_val;
  }
  

void SimInfo::read_simulation_params(const hid_t file_id)
  {
  using HDF5::read_string_set;

  // Open the data set

  if (file_id >= 0)
    {
    std::map<std::string, int> int_params;
    std::map<std::string, double> real_params;

    if (read_simulation_parameters(flash_version, file_id, int_params, 
				   real_params) < 0)
      throw Except("Unable to read simulations parameters", __FILE__, 
		   __LINE__);

    total_blocks = static_cast<unsigned int>(int_params[Total_Blocks_Name]);

    total_particles 
      = static_cast<unsigned int>(int_params[Total_Particles_Name]);

    sim_time = real_params[Time_Name];
    sim_dt = real_params[Timestep_Name];

    step_count = static_cast<unsigned int>(int_params[Step_Count_Name]);

    if (valid_mesh)
      {
      block_dims.resize(mesh_dims);

      switch (mesh_dims)
	{
	case 3 :
	  block_dims[2] = static_cast<unsigned int>(int_params[NZB_Name]);

	case 2 :
	  block_dims[1] = static_cast<unsigned int>(int_params[NYB_Name]);

	case 1 :
	  block_dims[0] = static_cast<unsigned int>(int_params[NXB_Name]);
	}

      // Count the number of cells per block

      block_cell_count = 1;

      for (unsigned int i = 0; i < mesh_dims; i++)
	block_cell_count *= block_dims[i];

      if (block_cell_count < 1)
	throw Except("No cells in block", __FILE__, __LINE__);
      }
    else
      {
      block_dims.clear();
      block_cell_count = 0;
      }
    }
  else
    throw Except("Invalid HDF5 file handle", __FILE__, __LINE__);
  }


void SimInfo::read_runtime_params(const hid_t file_id)
  {
  read_logical_runtime_params(file_id);
  read_integer_runtime_params(file_id);
  read_real_runtime_params(file_id);
  read_string_runtime_params(file_id);
  }


void SimInfo::read_logical_runtime_params(const hid_t file_id)
  {
  if (HDF5::read_param_table(file_id, Bool_Params_Name, 
			     logical_runtime_parameters) < 0)
    throw Except("Unable to read boolean runtime parameters", __FILE__,
		 __LINE__);
  }


void SimInfo::read_integer_runtime_params(const hid_t file_id)
  {
  if (HDF5::read_param_table(file_id, Int_Params_Name, 
			     integer_runtime_parameters) < 0)
    throw Except("Unable to read integer runtime parameters", __FILE__, 
		 __LINE__);

  if (valid_mesh)
    {
    base_block_dims.resize(mesh_dims);

    switch(mesh_dims)
      {
      case 3 :
	{
	const int base_blocks_z
	  = get_integer_runtime_param(BaseBlocks_Z_Name);
	  
	base_block_dims[2] = static_cast<unsigned int>(base_blocks_z);
	}
	// Continue to case 2

      case 2 :
	{
	const int base_blocks_y
	  = get_integer_runtime_param(BaseBlocks_Y_Name);

	base_block_dims[1] = static_cast<unsigned int>(base_blocks_y);
	}
	// Continue to case 1

      case 1 :
	{
	const int base_blocks_x
	  = get_integer_runtime_param(BaseBlocks_X_Name);

	base_block_dims[0] = static_cast<unsigned int>(base_blocks_x);
	}
	break;

      default :
	// ERROR CHECKING???
	break;
      }
    }
  else
    base_block_dims.clear();
  }


void SimInfo::read_real_runtime_params(const hid_t file_id)
  {
  if (HDF5::read_param_table(file_id, Float_Params_Name, 
			     real_runtime_parameters) < 0)
    throw Except("Unable to read real runtime parameters", __FILE__,
		 __LINE__);

  if (valid_mesh)
    {
    volume_min_bounds.resize(mesh_dims);
    volume_max_bounds.resize(mesh_dims);

    switch(mesh_dims)
      {
      case 3 :
	volume_min_bounds[2] = get_real_runtime_param(VolumeMin_Z_Name);
	volume_max_bounds[2] = get_real_runtime_param(VolumeMax_Z_Name);
	// Continue to case 2

      case 2 :
	volume_min_bounds[1] = get_real_runtime_param(VolumeMin_Y_Name);
	volume_max_bounds[1] = get_real_runtime_param(VolumeMax_Y_Name);
	// Continue to case 1

      case 1 :
	volume_min_bounds[0] = get_real_runtime_param(VolumeMin_X_Name);
	volume_max_bounds[0] = get_real_runtime_param(VolumeMax_X_Name);
	break;

      default :
	// ERROR CHECKING???
	break;
      }

    volume_center.resize(mesh_dims);
    volume_width.resize(mesh_dims);

    for (unsigned int i = 0; i < mesh_dims; i++)
      {
      const double axis_min_bound = volume_min_bounds[i];
      const double axis_max_bound = volume_max_bounds[i];

      volume_center[i] = 0.5 * (axis_min_bound + axis_max_bound);
      volume_width[i] = axis_max_bound - axis_min_bound;
      }
    }
  else
    {
    volume_min_bounds.clear();
    volume_max_bounds.clear();

    volume_center.clear();
    volume_width.clear();
    }
  }


void SimInfo::read_string_runtime_params(const hid_t file_id)
  {
  if (HDF5::read_param_table(file_id, String_Params_Name, 
			     string_runtime_parameters) < 0)
    throw Except("Unable to read string runtime parameters", __FILE__,
		 __LINE__);
  }


}  // End namespace File
}  // End namespace QuickFlash
