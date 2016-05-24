// C++ header file quickflash_file_dataset.hpp

/*
  By Nathan C. Hearn
     October 4, 2008

  Buffered access to HDF5 dataset.
*/


#ifndef QUICKFLASH_FILE_DATASET_HPP
#define QUICKFLASH_FILE_DATASET_HPP


#include <hdf5.h>
#include <time.h>
#include <vector>
#include <string>
#include "quickflash_block_blockdata.hpp"
#include "quickflash_cache_indexcache.hpp"


namespace QuickFlash
{
namespace File
{

class Dataset
  {
  public :

    typedef Block::BlockData<double> BlockData_Double;

  public :

    Dataset();

    Dataset(const hid_t file_container_id);

    ~Dataset();    

    bool ready() const { return dataset_ready; }

    const std::string & get_name() const { return name; }

    unsigned int get_num_blocks() const { return num_blocks; }
    unsigned int get_block_size() const { return block_size; }

    const std::vector<unsigned int> & get_block_dims() const 
      { return block_dims; }

    void reset();
    void reset(const hid_t file_container_id);

    void open_dataset(const std::string & dataset_name, 
		      const unsigned int space_dims,
		      const unsigned int buffer_length=1);

    void close_dataset();

    double get_min_value() const { return data_min_value; }
    double get_max_value() const { return data_max_value; }

    double get_file_read_time() const
      { 
      const double read_time_seconds 
	= static_cast<double>(accum_file_read_cpu_clock) 
	  * oneOver_clocks_per_sec;

      return read_time_seconds;
      }

    unsigned int get_buffer_size() const { return buffer_size; }

    unsigned int get_buffer_hits() const { return buffer_hits; }
    unsigned int get_buffer_misses() const { return buffer_misses; }

    void set_cache_size(const unsigned int max_cache_size=0,
			const unsigned int clean_cache_size=0) const
      {
      cache_size = max_cache_size;

      if (cache_size > 0)
	{
	if (cache.get_num_items() != num_blocks)
	  cache.reset(num_blocks);

	cache.set_max_cache_size(cache_size, clean_cache_size);
	}
      else
	cache.reset();
      }

    unsigned int get_cache_size() const { return cache_size; }

    unsigned int get_cache_stored_items() const 
      { return cache.get_num_cache_items(); }

    void clear_cache() const { cache.flush(); }

    unsigned int get_cache_hits() const { return cache_hits; }
    unsigned int get_cache_misses() const { return cache_misses; }

    void set_report_stats(const bool state=true) const 
      { report_stats = state; }

    void report_buffer_stats(bool report_if_unused=false) const;

    void report_cache_stats(bool report_if_unused=false) const;

    void report_all_stats(bool report_if_unused=false) const
      {
      report_buffer_stats(report_if_unused);
      report_cache_stats(report_if_unused);
      }

    void reset_buffer_stats() const;
    void reset_cache_stats() const;

    void reset_all_stats() const
      {
      reset_buffer_stats();
      reset_cache_stats();
      }

    void get_block_data(const unsigned int block_index, 
			std::vector<double> & data) const;

    void get_block_data(const unsigned int block_index,
			BlockData_Double & block) const
      {
      block.reset(block_dims);
      store_block_data(block_index, block);
      }

  private :

    void open_set(const std::string & dataset_name);
    void read_set_info(const unsigned int space_dims);
    void load_block(const unsigned int block_index) const;

    void store_block_data(const unsigned int block_index, 
			  BlockData_Double & block) const;

    const double * get_buffer_block_ptr(const unsigned int block_index) const;

  private :

    static const double oneOver_clocks_per_sec;

    hid_t container_id;  // Just a reference to an open HDF5 container

    hid_t dataset_id;

    bool dataset_ready;

    std::string name;

    unsigned int buffer_size;  // Number of blocks in buffer

    std::vector<unsigned int> dataset_dims;

    unsigned int block_size;   // Number of elements in block
    std::vector<unsigned int> block_dims;
    unsigned int num_blocks;  // Number of blocks in dataset

    double data_min_value;
    double data_max_value;

    mutable bool buffer_loaded;
    mutable unsigned int buffer_first_index;  // Index of first block in buffer
    mutable double * buffer;  // Stored in Fortran-style order
    mutable unsigned int buffer_hits;
    mutable unsigned int buffer_misses;

    mutable Cache::IndexCache< std::vector<double> > cache;  // C-style order
    mutable unsigned int cache_size;
    mutable unsigned int cache_hits;
    mutable unsigned int cache_misses;

    mutable bool report_stats;

    mutable clock_t accum_file_read_cpu_clock;
  };


}  // End namespace File
}  // End namespace QuickFlash


#endif  // QUICKFLASH_FILE_DATASET_HPP
