// C++ header file quickflash_hdf5.hpp

/*
  By Nathan C. Hearn
     September 15, 2006

  HDF5 utility functions
*/


#ifndef QUICKFLASH_HDF5_HPP
#define QUICKFLASH_HDF5_HPP


#include <hdf5.h>
#include <vector>
#include <string>
#include <map>


namespace QuickFlash
{
namespace HDF5
{

extern const unsigned int Parameter_Name_Size;
extern const unsigned int Parameter_String_Size;

extern const char Parameter_Name_Label[];
extern const char Parameter_Value_Label[];


// Functions

int get_dataset_dims(const hid_t dataset_id, std::vector<unsigned int> & dims);

int get_block_dataset_info(const hid_t dataset_id, unsigned int & num_blocks,
			   std::vector<unsigned int> & block_dims);

int get_block_dataset_info(const hid_t dataset_id, 
			   std::vector<unsigned int> & dataset_dims,
			   unsigned int & num_blocks,
			   std::vector<unsigned int> & block_dims);

int read_attribute(const hid_t container_id, const char attribute_name[],
		   double & value);

int read_attribute(const hid_t container_id, const char attribute_name[],
		   int & value);

int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector<double> & data);

int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector<int> & data);

int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector< std::vector<double> > & data);

int read_data(const hid_t container_id, const char dataset_name[],
	      const unsigned int num_vector_elems,
	      std::vector< std::vector<double> > & data, 
	      const double default_elem_value=0.0);

int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector< std::vector<int> > & data);

int read_data(const hid_t container_id, const char dataset_name[],
	      const unsigned int num_vector_elems,
	      std::vector< std::vector<int> > & data,
	      const int default_elem_value=0);

int read_data_stripe(const hid_t container_id, const char dataset_name[],
		     const unsigned int column_axis, 
		     const std::vector<unsigned int> & column_coords,
		     std::vector<double> & data);

void strip_trailing_spaces(const unsigned int str_len, char str[]);
void strip_trailing_spaces(std::string & str);

int read_param_table(const hid_t container_id, const char table_name[],
		     std::map<std::string, bool> & params);

int read_param_table(const hid_t container_id, const char table_name[],
		     std::map<std::string, int> & params);

int read_param_table(const hid_t container_id, const char table_name[],
		     std::map<std::string, double> & params);

int read_param_table(const hid_t container_id, const char table_name[],
		     std::map<std::string, std::string> & params);

int read_compound_param(const hid_t container_id, const char table_name[],
			const char elem_name[], bool & value);

int read_compound_param(const hid_t container_id, const char table_name[],
			const char elem_name[], int & value);

int read_compound_param(const hid_t container_id, const char table_name[],
			const char elem_name[], double & value);

// ADD read_compound_param FOR STRINGS -- BE CAREFUL ABOUT LENGTHS!!!

int read_string_set(const hid_t container_id, const char dataset_name[],
		    std::vector<std::string> & data);

int get_num_members(const hid_t container_id);

int member_index(const hid_t container_id, const std::string & name);
int member_index(const hid_t container_id, const char name[]);

bool member_present(const hid_t container_id, const std::string & name);
bool member_present(const hid_t container_id, const char name[]);


}  // End namespace HDF5
}  // End namespace QuickFlash


#endif  // QUICKFLASH_HDF5_HPP
