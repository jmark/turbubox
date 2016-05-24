// C++ program file quickflash_file_flashdefs.cpp

/*
  By Nathan C. Hearn
     December 13, 2006

  Flash datafile definitions.  Flash3 detection as suggested by Dean Townsley.
*/


#include "quickflash_file_flashdefs.hpp"
#include <hdf5.h>
#include <string>
#include <vector>
#include <map>
#include "quickflash_hdf5.hpp"
#include "quickflash_except.hpp"


namespace QuickFlash
{
namespace File
{

const unsigned int Flash2_FileVersion_Min = 7;
const unsigned int Flash2_FileVersion_Max = 7;

const unsigned int Flash3_FileVersion_Min = 8;
const unsigned int Flash3_FileVersion_Max = 9;

const char File_Format_Name[] = "file format version";

const char Sim_Params_Name_Flash2[] = "simulation parameters";
const char Sim_Info_Name_Flash3[] = "sim info";

const char Int_Params_Name[] = "integer runtime parameters";
const char Bool_Params_Name[] = "logical runtime parameters";
const char Float_Params_Name[] = "real runtime parameters";
const char String_Params_Name[] = "string runtime parameters";

const char Int_Scalars_Name_Flash3[] = "integer scalars";
const char Float_Scalars_Name_Flash3[] = "real scalars";

const char Dimensionality_Name_Flash3[] = "dimensionality";

const unsigned int Dimensionality_Min = 1;
const unsigned int Dimensionality_Max = 3;

const char Total_Blocks_Name[] = "total_blocks";
const char Total_Blocks_Name_Flash2[] = "total blocks";
const char Total_Blocks_Name_Flash3[] = "globalnumblocks";

const char Total_Particles_Name[] = "total_particles";
const char Total_Particles_Name_Flash3[] = "globalnumparticles";

const char Time_Name[] = "time";

const char Timestep_Name[] = "dt";
const char Timestep_Name_Flash2[] = "timestep";
const char Timestep_Name_Flash3[] = "dt";

const char Step_Count_Name[] = "timestep_index";
const char Step_Count_Name_Flash2[] = "number of steps";
const char Step_Count_Name_Flash3[] = "nstep";

const char NXB_Name[] = "nxb";
const char NYB_Name[] = "nyb";
const char NZB_Name[] = "nzb";

const char BaseBlocks_X_Name[] = "nblockx";
const char BaseBlocks_Y_Name[] = "nblocky";
const char BaseBlocks_Z_Name[] = "nblockz";

const char VolumeMin_X_Name[] = "xmin";
const char VolumeMin_Y_Name[] = "ymin";
const char VolumeMin_Z_Name[] = "zmin";
const char VolumeMax_X_Name[] = "xmax";
const char VolumeMax_Y_Name[] = "ymax";
const char VolumeMax_Z_Name[] = "zmax";

const char Geometry_Label[] = "geometry";

const char Tree_Struct_Name[] = "gid";

const unsigned int Tree_Struct_Elems_1d = 5;
const unsigned int Tree_Struct_Elems_2d = 9;
const unsigned int Tree_Struct_Elems_3d = 15;

const char Coordinates_Name[] = "coordinates";
const char BlockWidth_Name[] = "block size";
const char RefineLevel_Name[] = "refine level";
const char BlockType_Name[] = "node type";
const char ProcessID_Name[] = "processor number";

const char Variables_Name[] = "unknown names";

const unsigned int Variable_Name_Length_v7 = 4;
const unsigned int Variable_Name_Length_v8 = 4;
const unsigned int Variable_Name_Length_v9 = 4;

const char Variable_Name_Pad_Char = ' ';

const char Particle_Variables_Name_Flash3[] = "particle names";

const char Particle_Variables_Tag_Name[] = "tag";

const char Particle_Data_Name[] = "tracer particles";

const char Dataset_Attribute_MinValue_Name[] = "minimum";
const char Dataset_Attribute_MaxValue_Name[] = "maximum";


int get_flash_version(const hid_t file_id, unsigned int & file_version,
		      FlashVersion & flash_version)
  {
  int ret_val = -1;  // Error state

  // See if we are dealing with Flash 2 or 3

  if (HDF5::member_present(file_id, File_Format_Name))
    {
    flash_version = Flash2;

    std::vector<int> file_format_data;

    if (HDF5::read_data(file_id, File_Format_Name, file_format_data) >= 0)
      {
      if (file_format_data.size() == 1)
	{
	file_version = static_cast<unsigned int>(file_format_data[0]);

	if ((file_version >= Flash2_FileVersion_Min)
	    && (file_version <= Flash2_FileVersion_Max))
	  ret_val = 0;  // OK
	}
      }
    }
  else if (HDF5::member_present(file_id, Sim_Info_Name_Flash3))
    {
    flash_version = Flash3;

    int raw_version_number;

    if (HDF5::read_compound_param(file_id, Sim_Info_Name_Flash3, 
				  File_Format_Name, raw_version_number) >= 0)
      {
      file_version = raw_version_number;

      if ((file_version >= Flash3_FileVersion_Min)
	  && (file_version <= Flash3_FileVersion_Max))
	ret_val = 0;  // OK
      }
    }

  return ret_val;
  }


int get_variable_name_length(const unsigned int file_version, 
			     unsigned int & variable_name_length)
  {
  int ret_val = 0;  // OK

  switch (file_version)
    {
    case 7 :
      variable_name_length = Variable_Name_Length_v7;
      break;

    case 8 :
      variable_name_length = Variable_Name_Length_v8;
      break;

    case 9 :
      variable_name_length = Variable_Name_Length_v9;
      break;

    default :
      ret_val = -1;  // Error state
      break;
    }

  return ret_val;
  }


int pad_variable_name(const unsigned int variable_name_length, 
		      const std::string & orig_name, std::string & new_name)
  { 
  int ret_val = -1;  // Error state (may mean empty orig_name)

  const unsigned int str_len = orig_name.size();

  new_name = orig_name;

  new_name.resize(variable_name_length, Variable_Name_Pad_Char);

  if (str_len > variable_name_length)
    ret_val = 1;  // Original name exceeds length
  else if (str_len > 0)
    ret_val = 0;  // OK    

  return ret_val;
  }


int pad_variable_name_file_version(const unsigned int file_version,
				   const std::string & orig_name,
				   std::string & new_name)
  {
  int ret_val = -1;  // Error state

  unsigned int var_name_len = 0;

  if (get_variable_name_length(file_version, var_name_len) >= 0)
    ret_val = pad_variable_name(var_name_len, orig_name, new_name);

  return ret_val;
  }


int get_dimensionality_tree_struct(const unsigned int tree_struct_elems,
				   unsigned int & dims)
  {
  int ret_val = -1;  // Error state

  if (tree_struct_elems == Tree_Struct_Elems_1d)
    {
    dims = 1;
    ret_val = 0;  // OK
    }
  else if (tree_struct_elems == Tree_Struct_Elems_2d)
    {
    dims = 2;
    ret_val = 0;  // OK
    }
  else if (tree_struct_elems == Tree_Struct_Elems_3d)
    {
    dims = 3;
    ret_val = 0;  // OK
    }

  return ret_val;
  }


int get_mesh_dimensionality(const hid_t file_id, 
			    const unsigned int file_version, 
			    unsigned int & dims)
  {
  int ret_val = -1;  // Error state

  switch (file_version)
    {
    case 7 :
    case 8 :

      if (HDF5::member_present(file_id, Tree_Struct_Name))
	{
	std::vector<unsigned int> gid_dims;

	hid_t dataset_id = H5Dopen(file_id, Tree_Struct_Name);

	if (dataset_id >= 0)
	  {
	  if (HDF5::get_dataset_dims(dataset_id, gid_dims) >= 0)
	    {
	    if (gid_dims.size() == 2)
	      {
	      const unsigned int tree_struct_elems = gid_dims[1];

	      ret_val
		= get_dimensionality_tree_struct(tree_struct_elems, dims);
	      }
	    }

	  H5Dclose(dataset_id);
	  dataset_id = 0;
	  }
	}

      break;

    case 9 :

      if (HDF5::member_present(file_id, Tree_Struct_Name))
	{
	std::map<std::string, int> int_scalars;

	if (HDF5::read_param_table(file_id, Int_Scalars_Name_Flash3, 
				   int_scalars) >= 0)
	  {
	  const std::map<std::string, int>::const_iterator dims_iter
	    = int_scalars.find(std::string(Dimensionality_Name_Flash3));

	  if (dims_iter != int_scalars.end())
	    {
	    dims = static_cast<unsigned int>(dims_iter->second);

	    if ((dims >= Dimensionality_Min) && (dims <= Dimensionality_Max))
	      ret_val = 0;  // OK
	    }
	  else
	    {
	    // Default to the version 8 way

	    std::vector<unsigned int> gid_dims;

	    hid_t dataset_id = H5Dopen(file_id, Tree_Struct_Name);

	    if (dataset_id >= 0)
	      {
	      if (HDF5::get_dataset_dims(dataset_id, gid_dims) >= 0)
		{
		if (gid_dims.size() == 2)
		  {
		  const unsigned int tree_struct_elems = gid_dims[1];

		  ret_val = get_dimensionality_tree_struct(tree_struct_elems, 
							   dims);
		  }
		}

	      H5Dclose(dataset_id);
	      dataset_id = 0;
	      }
	    }
	  }
	}

      break;

    default :

      break;  // Retain error state
    }

  return ret_val;
  }


int read_simulation_parameters(const FlashVersion flash_code_version,
			       const hid_t file_id,
			       std::map<std::string, int> & int_param_map,
			       std::map<std::string, double> & real_param_map)
  {
  int ret_val = -1;  // Error state

  switch (flash_code_version)
    {
    case Flash2 :
      {
      hid_t dataset_id = H5Dopen(file_id, Sim_Params_Name_Flash2);

      if (dataset_id >= 0)
	{
	ret_val = read_int_simulation_params_flash2(dataset_id, int_param_map);

	if (ret_val >= 0)
	  {
	  ret_val = read_real_simulation_params_flash2(dataset_id, 
						       real_param_map);

	  if (ret_val >= 0)
	    ret_val = 0;  // OK
	  }

	H5Dclose(dataset_id);
	dataset_id = -1;
	}
      }
      break;

    case Flash3 :
      {
      ret_val = read_int_simulation_params_flash3(file_id, int_param_map);

      if (ret_val >= 0)
	{
	ret_val = read_real_simulation_params_flash3(file_id, real_param_map);

	if (ret_val >= 0)
	  ret_val = 0;  // OK
	}
      }
      break;

    default :
      // Maintain error state
      break;
    }

  return ret_val;
  }


int read_int_simulation_params_flash2(const hid_t dataset_id,
				      std::map<std::string, int> & param_map)
  {
  int ret_val = -1;  // Error state

  // Read the integer components

  const int num_ints = 5;
  const size_t num_int_bytes = sizeof(int) * num_ints;

  hid_t datatype_int_id = H5Tcreate(H5T_COMPOUND, num_int_bytes);

  if (datatype_int_id >= 0)
    {
    int i = 0;

    H5Tinsert(datatype_int_id, Total_Blocks_Name_Flash2, sizeof(int) * i, 
	      H5T_NATIVE_INT);
    i++;

    H5Tinsert(datatype_int_id, Step_Count_Name_Flash2, sizeof(int) * i, 
	      H5T_NATIVE_INT);
    i++;

    H5Tinsert(datatype_int_id, NXB_Name, sizeof(int) * i, H5T_NATIVE_INT);
    i++;

    H5Tinsert(datatype_int_id, NYB_Name, sizeof(int) * i, H5T_NATIVE_INT);
    i++;

    H5Tinsert(datatype_int_id, NZB_Name, sizeof(int) * i, H5T_NATIVE_INT);
    i++;

    int int_array[num_ints];

    herr_t status = H5Dread(dataset_id, datatype_int_id, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, int_array);

    if (status >= 0)
      {
      i = 0;

      param_map[Total_Blocks_Name] = int_array[i++];
      param_map[Step_Count_Name] = int_array[i++];

      param_map[NXB_Name] = int_array[i++];
      param_map[NYB_Name] = int_array[i++];
      param_map[NZB_Name] = int_array[i++];

      ret_val = 0;  // OK
      }

    // Particles -- CURRENTLY UNSUPPPORTED FOR FLASH2!!!!!

    param_map[Total_Particles_Name] = 0;  // For now ...

    H5Tclose(datatype_int_id);
    datatype_int_id = -1;
    }

  return ret_val;
  }


int read_real_simulation_params_flash2(const hid_t dataset_id,
				       std::map<std::string, double> 
				         & param_map)
  {
  int ret_val = -1;  // Error state

  // Read the floating point components

  const int num_floats = 2;
  const size_t num_float_bytes = sizeof(float) * num_floats;

  hid_t datatype_float_id = H5Tcreate(H5T_COMPOUND, num_float_bytes);

  if (datatype_float_id >= 0)
    {
    int i = 0;

    H5Tinsert(datatype_float_id, Time_Name, sizeof(float) * i, 
              H5T_NATIVE_FLOAT);
    i++;

    H5Tinsert(datatype_float_id, Timestep_Name_Flash2, sizeof(float) * i, 
	      H5T_NATIVE_FLOAT);
    i++;

    float float_array[num_floats];

    herr_t status = H5Dread(dataset_id, datatype_float_id, H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, float_array);

    if (status >= 0)
      {
      i = 0;

      param_map[Time_Name] = static_cast<double>(float_array[i++]);
      param_map[Timestep_Name] = static_cast<double>(float_array[i++]);

      ret_val = 0;  // OK
      }

    H5Tclose(datatype_float_id);
    datatype_float_id = -1;
    }

  return ret_val;
  }


int read_int_simulation_params_flash3(const hid_t file_id,
				      std::map<std::string, int> & param_map)
  {
  int ret_val = -1;  // Error state

  // NOTE: dimensions must be properly defined before calling this function

  std::map<std::string, int> int_params;

  if (HDF5::read_param_table(file_id, Int_Scalars_Name_Flash3, int_params) 
      >= 0)
    {
    // WE SHOULD CHECK PRESENCE OF PARAMETERS IN MAP!!!

    param_map[Total_Blocks_Name] = int_params[Total_Blocks_Name_Flash3];
    param_map[Step_Count_Name] = int_params[Step_Count_Name_Flash3];

    param_map[NXB_Name] = int_params[NXB_Name];
    param_map[NYB_Name] = int_params[NYB_Name];
    param_map[NZB_Name] = int_params[NZB_Name];

    const std::map<std::string, int>::const_iterator part_num_ptr 
      = int_params.find(Total_Particles_Name_Flash3);

    int total_particles = 0;

    if (part_num_ptr != int_params.end())
      total_particles = part_num_ptr->second;

    param_map[Total_Particles_Name] = total_particles;

    ret_val = 0;  // OK
    }

  return ret_val;
  }


int read_real_simulation_params_flash3(const hid_t file_id,
				       std::map<std::string, double> 
				         & param_map)
  {
  int ret_val = -1;

  std::map<std::string, double> float_params;

  if (HDF5::read_param_table(file_id, Float_Scalars_Name_Flash3, float_params) 
      >= 0)
    {
    // WE SHOULD CHECK PRESENCE OF PARAMETERS IN MAP!!!

    param_map[Time_Name] = float_params[Time_Name];
    param_map[Timestep_Name] = float_params[Timestep_Name_Flash3];

    ret_val = 0;
    }

  return ret_val;
  }


}  // End namespace File
}  // End namespace QuickFlash
