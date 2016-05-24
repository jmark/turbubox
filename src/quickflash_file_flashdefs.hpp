// C++ header file quickflash_file_flashdefs.hpp

/*
  By Nathan C. Hearn
     December 13, 2006

  Flash datafile definitions.
*/


#ifndef QUICKFLASH_FILE_FLASHDEFS_HPP
#define QUICKFLASH_FILE_FLASHDEFS_HPP


#include <hdf5.h>
#include <string>
#include <map>


namespace QuickFlash
{
namespace File
{

enum FlashVersion { FlashUnknown, Flash2, Flash3 };

extern const unsigned int Flash2_FileVersion_Min;
extern const unsigned int Flash2_FileVersion_Max;

extern const unsigned int Flash3_FileVersion_Min;
extern const unsigned int Flash3_FileVersion_Max;

extern const char File_Format_Name[];

extern const char Sim_Params_Name_Flash2[];
extern const char Sim_Info_Name_Flash3[];

extern const char Int_Params_Name[];
extern const char Bool_Params_Name[];
extern const char Float_Params_Name[];
extern const char String_Params_Name[];

extern const char Int_Scalars_Name_Flash3[];
extern const char Float_Scalars_Name_Flash3[];

extern const char Dimensionality_Name_Flash3[];

extern const unsigned int Dimensionality_Min;
extern const unsigned int Dimensionality_Max;

extern const char Total_Blocks_Name[];
extern const char Total_Blocks_Name_Flash2[];
extern const char Total_Blocks_Name_Flash3[];

extern const char Total_Particles_Name[];
extern const char Total_Particles_Name_Flash3[];

extern const char Time_Name[];

extern const char Timestep_Name[];
extern const char Timestep_Name_Flash2[];
extern const char Timestep_Name_Flash3[];

extern const char Step_Count_Name[];
extern const char Step_Count_Name_Flash2[];
extern const char Step_Count_Name_Flash3[];

extern const char NXB_Name[];
extern const char NYB_Name[];
extern const char NZB_Name[];

extern const char BaseBlocks_X_Name[];
extern const char BaseBlocks_Y_Name[];
extern const char BaseBlocks_Z_Name[];

extern const char VolumeMin_X_Name[];
extern const char VolumeMin_Y_Name[];
extern const char VolumeMin_Z_Name[];
extern const char VolumeMax_X_Name[];
extern const char VolumeMax_Y_Name[];
extern const char VolumeMax_Z_Name[];

extern const char Geometry_Label[];

extern const char Tree_Struct_Name[];

extern const unsigned int Tree_Struct_Elems_1d;
extern const unsigned int Tree_Struct_Elems_2d;
extern const unsigned int Tree_Struct_Elems_3d;

extern const char Coordinates_Name[];
extern const char BlockWidth_Name[];
extern const char RefineLevel_Name[];
extern const char BlockType_Name[];
extern const char ProcessID_Name[];

extern const char Variables_Name[];

extern const unsigned int Variable_Name_Length_v7;
extern const unsigned int Variable_Name_Length_v8;
extern const unsigned int Variable_Name_Length_v9;

extern const char Variable_Name_Pad_Char;

extern const char Particle_Variables_Name_Flash3[];

extern const char Particle_Variables_Tag_Name[];

extern const char Particle_Data_Name[];

extern const char Dataset_Attribute_MinValue_Name[];
extern const char Dataset_Attribute_MaxValue_Name[];


int get_flash_version(const hid_t file_id, unsigned int & file_version,
		      FlashVersion & flash_version);

int get_variable_name_length(const unsigned int file_version, 
			     unsigned int & variable_name_length);

int pad_variable_name(const unsigned int variable_name_length,
		      const std::string & orig_name, std::string & new_name);

int pad_variable_name_file_version(const unsigned int file_version,
				   const std::string & orig_name,
				   std::string & new_name);

int get_dimensionality_tree_struct(const unsigned int tree_struct_elems,
				   unsigned int & dims);

int get_mesh_dimensionality(const hid_t file_id, 
			    const unsigned int file_version, 
			    unsigned int & dims);

int read_simulation_parameters(const FlashVersion flash_code_version,
			       const hid_t file_id,
			       std::map<std::string, int> & int_param_map,
			       std::map<std::string, double> & real_param_map);

int read_int_simulation_params_flash2(const hid_t dataset_id,
				      std::map<std::string, int> & param_map);

int read_real_simulation_params_flash2(const hid_t dataset_id,
				       std::map<std::string, double> 
				         & param_map);

int read_int_simulation_params_flash3(const hid_t file_id,
				      std::map<std::string, int> & param_map);


int read_real_simulation_params_flash3(const hid_t file_id,
				       std::map<std::string, double> 
				         & param_map);


}  // End namespace File
}  // End namespace QuickFlash


#endif  // QUICKFLASH_FILE_FLASHDEFS_HPP
