// C++ header file quickflash_file_siminfo.hpp

/*
  By Nathan C. Hearn
     October 4, 2008

  Container for general simulation information.
*/


#ifndef QUICKFLASH_FILE_SIMINFO_HPP
#define QUICKFLASH_FILE_SIMINFO_HPP


#include <hdf5.h>
#include <vector>
#include <string>
#include <map>
#include "quickflash_file_flashdefs.hpp"
#include "quickflash_geometry.hpp"


namespace QuickFlash
{
namespace File
{

class SimInfo
  {
  public :

    SimInfo();

    SimInfo(const hid_t file_id, const FlashVersion flash_code_version,
	    const unsigned int flash_file_version, const bool skip_mesh=false);

    SimInfo(const SimInfo & source);

    ~SimInfo() { }

    SimInfo & operator=(const SimInfo & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const hid_t file_id, const FlashVersion flash_code_version,
	       const unsigned int flash_file_version, 
	       const bool skip_mesh=false);

    void reset(const SimInfo & source);

    FlashVersion get_flash_code_version() const { return flash_version; }
    unsigned int get_flash_file_version() const { return file_version; }

    bool mesh_present() const { return valid_mesh; }
    bool particles_present() const { return valid_particles; }

    unsigned int get_mesh_dims() const { return mesh_dims; }

    Geometry::GeometryType get_mesh_geometry() const { return geom_type; }

    unsigned int get_num_blocks() const { return total_blocks; }
    unsigned int get_num_particles() const { return total_particles; }

    const std::vector<unsigned int> & get_base_block_dims() const
      { return base_block_dims; }

    const std::vector<double> & get_volume_minbounds() const
      { return volume_min_bounds; }

    const std::vector<double> & get_volume_maxbounds() const
      { return volume_max_bounds; }

    const std::vector<double> & get_volume_center() const
      { return volume_center; }

    const std::vector<double> & get_volume_width() const
      { return volume_width; }

    double get_sim_time() const { return sim_time; }
    double get_sim_dt() const { return sim_dt; }

    unsigned int get_step_count() const { return step_count; }

    const std::vector<unsigned int> & get_block_dims() const 
      { return block_dims; }

    unsigned int get_cells_per_block() const { return block_cell_count; }

    bool logical_runtime_param_present(const std::string & name) const;
    bool get_logical_runtime_param(const std::string & name) const;

    bool integer_runtime_param_present(const std::string & name) const;
    int get_integer_runtime_param(const std::string & name) const;

    bool real_runtime_param_present(const std::string & name) const;
    double get_real_runtime_param(const std::string & name) const;

    bool string_runtime_param_present(const std::string & name) const;
    const std::string & get_string_runtime_param(const std::string & name) 
      const;

    unsigned int get_num_mesh_variables() const 
      { return num_mesh_variables; }

    const std::vector<std::string> & get_mesh_variable_names() const 
      { return mesh_variable_names; }

    bool mesh_variable_present(const std::string & variable_name) const
      { return mesh_variable_present(variable_name.c_str()); }

    bool mesh_variable_present(const char variable_name[]) const;

  private :

    void read_simulation_params(const hid_t file_id);
    
    void read_runtime_params(const hid_t file_id);
    void read_logical_runtime_params(const hid_t file_id);
    void read_integer_runtime_params(const hid_t file_id);
    void read_real_runtime_params(const hid_t file_id);
    void read_string_runtime_params(const hid_t file_id);

  private :

    FlashVersion flash_version;
    unsigned int file_version;

    bool valid_mesh;
    bool valid_particles;

    unsigned int mesh_dims;

    unsigned int total_blocks;
    unsigned int total_particles;

    double sim_time;
    double sim_dt;

    unsigned int step_count;

    std::vector<unsigned int> block_dims;
    unsigned int block_cell_count;

    unsigned int num_mesh_variables;
    std::vector<std::string> mesh_variable_names;

    std::vector<unsigned int> base_block_dims;

    std::vector<double> volume_min_bounds;
    std::vector<double> volume_max_bounds;

    std::vector<double> volume_center;
    std::vector<double> volume_width;

    std::map<std::string, bool> logical_runtime_parameters;
    std::map<std::string, int> integer_runtime_parameters;
    std::map<std::string, double> real_runtime_parameters;
    std::map<std::string, std::string> string_runtime_parameters;

    Geometry::GeometryType geom_type;
  };


}  // End namespace File
}  // End namespace QuickFlash


#endif  // QUICKFLASH_FILE_SIMINFO_HPP
