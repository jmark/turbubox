// C++ program file quickflash_dataset.cpp

/*
  By Nathan C. Hearn
     September 15, 2006

  Buffered access to HDF5 dataset.

  NOTE: Buffered data is stored as in file, with Fortran ordering.  Cached
  data (and any direct storage into std::vector or BlockData types) is stored
  in normal C-style ordering.
*/


#include "quickflash_file_dataset.hpp"
#include <hdf5.h>
#include <time.h>
#include <string>
#include <iostream>
#include "quickflash_except.hpp"
#include "quickflash_hdf5.hpp"
#include "quickflash_counters.hpp"
#include "quickflash_file_flashdefs.hpp"

// Class Dataset


namespace QuickFlash
{
namespace File
{

// Class Dataset

const double Dataset::oneOver_clocks_per_sec 
  = 1.0 / static_cast<double>(CLOCKS_PER_SEC);  // Hope it's set by compiler


Dataset::Dataset() : 
  container_id(-1), dataset_id(-1), dataset_ready(false), name(), 
  buffer_size(0), dataset_dims(), block_size(0), block_dims(), 
  num_blocks(0), data_min_value(0.0), data_max_value(0.0),
  buffer_loaded(false), buffer_first_index(0), buffer(0), buffer_hits(0),
  buffer_misses(0), cache(), cache_size(0), cache_hits(0), cache_misses(0),
  report_stats(false), accum_file_read_cpu_clock(0)
  { }


Dataset::Dataset(const hid_t file_container_id) :
  container_id(-1), dataset_id(-1), dataset_ready(false), name(), 
  buffer_size(0), dataset_dims(), block_size(0), block_dims(), 
  num_blocks(0), data_min_value(0.0), data_max_value(0.0),
  buffer_loaded(false), buffer_first_index(0), buffer(0), buffer_hits(0),
  buffer_misses(0), cache(), cache_size(0), cache_hits(0), cache_misses(0),
  report_stats(false), accum_file_read_cpu_clock(0)
  { reset(file_container_id); }


Dataset::~Dataset() 
  { reset(); }


void Dataset::reset()
  {
  close_dataset();

  container_id = -1;  // Don't delete the object!  Someone else has it!
  }


void Dataset::reset(const hid_t file_container_id)
  {
  reset();

  container_id = file_container_id;

  // Reset name???
  }


void Dataset::open_dataset(const std::string & dataset_name, 
			   const unsigned int space_dims, 
			   const unsigned int buffer_length)
  {
  open_set(dataset_name);

  read_set_info(space_dims);

  // Determine buffer size

  unsigned int real_buffer_length = buffer_length;

  if (real_buffer_length < 1)
    real_buffer_length = 1;
  else if (real_buffer_length > num_blocks)
    real_buffer_length = num_blocks;

  const unsigned int num_elems = real_buffer_length * block_size;

  buffer = new double[num_elems];

  if (buffer == 0)
    throw Except("Buffer could not be allocated", __FILE__, __LINE__);

  buffer_size = real_buffer_length;

  cache.reset(num_blocks, 0);
  cache_size = 0;

  buffer_hits = 0;
  buffer_misses = 0;

  cache_hits = 0;
  cache_misses = 0;

  dataset_ready = true;
  }


void Dataset::close_dataset()
  {
  if (report_stats)
    report_all_stats();

  if (buffer != 0)
    {
    delete[] buffer;
    buffer = 0;
    }

  if (dataset_id >= 0)
    {
    H5Dclose(dataset_id);
    dataset_id = -1;
    }

  name = "";

  dataset_ready = false;

  buffer_size = 0;

  buffer_loaded = false;

  block_size = 0;
  block_dims.clear();

  num_blocks = 0;

  data_min_value = 0.0;
  data_max_value = 0.0;

  cache.reset();
  cache_size = 0;

  accum_file_read_cpu_clock = static_cast<clock_t>(0);
  }


void Dataset::report_buffer_stats(bool report_if_unused) const
  {
  using std::cerr;
  using std::endl;

  bool print_stats = false;

  if (report_if_unused)
    print_stats = true;
  else if ((buffer_hits > 0) || (buffer_misses > 0))
    print_stats = true;

  if (print_stats)
    {
    const double read_time_seconds = get_file_read_time();

    cerr << "Dataset " << name << " buffer hits [ " << buffer_hits
	 << " ] misses [ " << buffer_misses << " ] file read time [ " 
	 << read_time_seconds << " ]" << endl;
    }
  }


void Dataset::report_cache_stats(bool report_if_unused) const
  {
  using std::cerr;
  using std::endl;

  bool print_stats = false;

  if (report_if_unused)
    print_stats = true;
  else if ((cache_hits > 0) || (cache_misses > 0))
    print_stats = true;

  if (print_stats)
    {
    const unsigned int clean_count = cache.get_clean_counter();

    const unsigned int create_count = cache.get_num_created_items();

    const unsigned int multi_create_count 
      = cache.get_num_multi_created_items();

    cerr << "Dataset " << name << " cache hits [ " << cache_hits
	 << " ] misses [ " << cache_misses << " ] cleans [ " 
	 << clean_count << " ] creations [ " << create_count 
	 << " ] multi-creations [ " << multi_create_count << " ]" 
	 << endl;
    }
  }


void Dataset::get_block_data (
    const unsigned int block_index, 
	std::vector<double> & data
) const
{
  // If data is in the cache, use it; otherwise get the data from the buffer

    if (cache.data_present(block_index))
    {
        data = cache[block_index];
        cache_hits++;
        //std::cerr << "block_index = " << block_index << " (using cache)" << std::endl;
    }
    else
    {
        //std::cerr << "block_index = " << block_index << " (not using cache)" << std::endl;
        // Get a pointer to the block

        const double * block_ptr = get_buffer_block_ptr(block_index);

    // Copy the data -- correct for Fortran ordering

    data.resize(block_size);

    VectorCounter counter(block_dims);

    while (counter.in_bounds())
      {
      const unsigned int fortran_index = counter.get_scalar_index_reverse();
      const unsigned int local_index = counter.get_scalar_index();

      data[local_index] = block_ptr[fortran_index];

      counter.increment();
      }

    // If the cache is used, add the data to the cache

    if (cache_size > 0)
      {
      cache.add_data(block_index, data);
      cache_misses++;
      }
    }
  }

void Dataset::reset_buffer_stats() const
  {
  buffer_hits = 0;
  buffer_misses = 0;
  }


void Dataset::reset_cache_stats() const
  {
  cache_hits = 0;
  cache_misses = 0;

  cache.reset_access_counters();
  cache.reset_create_counters();
  }


void Dataset::open_set(const std::string & dataset_name)
  {
  close_dataset();

  if (container_id < 0)
    throw Except("Invalid container",  __FILE__, __LINE__);

  dataset_id = H5Dopen(container_id, dataset_name.c_str());

  if (dataset_id < 0)
    throw Except("Unable to open dataset", __FILE__, __LINE__);

  name = dataset_name;
  }


void Dataset::read_set_info(const unsigned int space_dims)
  {
  using HDF5::get_block_dataset_info;
  using HDF5::read_attribute;

  if (dataset_id < 0)
    throw Except("Dataset not open", __FILE__, __LINE__);

  // Read the dimensions from the dataspace

  if (get_block_dataset_info(dataset_id, dataset_dims, num_blocks, 
			     block_dims) < 0)
    throw Except("Unable to read dataset dimensions", __FILE__, __LINE__);

  // Reset block_dims to the dimensions of the simulation volume

  const unsigned int num_raw_block_dims = block_dims.size();

  if (num_raw_block_dims < 1)
    throw Except("Blocks have zero size", __FILE__, __LINE__);

  for (unsigned int i = space_dims; i < num_raw_block_dims; i++)
    if (block_dims[i] > 1)
      throw Except("File block dimensions incompatible with space dimensions",
		   __FILE__, __LINE__);

  block_dims.resize(space_dims, 1);

  // Count the number of cells in the block

  block_size = 1;

  for (unsigned int i = 0; i < space_dims; i++)
    block_size *= block_dims[i];

  // Get the minimum and maximum bounds

  if (read_attribute(dataset_id, Dataset_Attribute_MinValue_Name, 
		     data_min_value) < 0)
    throw Except("Unable to read minimum data value", __FILE__, __LINE__);

  if (read_attribute(dataset_id, Dataset_Attribute_MaxValue_Name, 
		     data_max_value) < 0)
    throw Except("Unable to read maximum data value", __FILE__, __LINE__);
  }


void Dataset::load_block(const unsigned int block_index) const
  {
  if (!dataset_ready)
    throw Except("Dataset is not available", __FILE__, __LINE__);

  if (block_index >= num_blocks)
    throw Except("Block index out of range", __FILE__, __LINE__);

  // Determine if blocks need to be loaded

  bool load_required = true;

  if (buffer_loaded)
    if ((block_index >= buffer_first_index) 
	&& (block_index < (buffer_first_index + buffer_size)))
      load_required = false;

  if (!load_required)
    buffer_hits++;
  else
    {
    buffer_misses++;

    bool error = false;

    const unsigned int max_first_buffer_index = num_blocks - buffer_size;

    unsigned int new_buffer_first_index = block_index;

    if (new_buffer_first_index > max_first_buffer_index)
      new_buffer_first_index = max_first_buffer_index;

    // Construct the data space for reading

    const clock_t begin_cpu_clock = clock();  // Time at start

    hid_t file_space_id = H5Dget_space(dataset_id);

    if (file_space_id < 0)
      error = true;
    else
      {
      const unsigned int space_dims = block_dims.size();
      const unsigned int hdf5_dims = dataset_dims.size();

      if (hdf5_dims < (space_dims + 1))
	error = true;
      else
	{
	hsize_t * start_array = new hsize_t[hdf5_dims];
	hsize_t * stride_array = new hsize_t[hdf5_dims];
	hsize_t * count_array = new hsize_t[hdf5_dims];
	hsize_t * block_array = new hsize_t[hdf5_dims];

	start_array[0] = static_cast<hsize_t>(new_buffer_first_index);
	stride_array[0] = 1;
	count_array[0] = static_cast<hsize_t>(buffer_size);
	block_array[0] = 1;

	const unsigned int space_dims_start = hdf5_dims - space_dims;

	for (unsigned int i = 1; i < hdf5_dims; i++)
	  {
	  unsigned int dim_length = 1;

	  if (i >= space_dims_start)
	    {
	    // Adjust for Fortran ordering...

	    const unsigned space_dims_index = hdf5_dims - (i + 1);

	    dim_length = block_dims[space_dims_index];
	    }

	  start_array[i] = 0;
	  stride_array[i] = 1;
	  count_array[i] = static_cast<hsize_t>(dim_length);
	  block_array[i] = 1;
	  }

	if (H5Sselect_hyperslab(file_space_id, H5S_SELECT_AND, start_array,
				stride_array, count_array, block_array) < 0)
	  error = true;
	else
	  {
	  hid_t mem_space_id = H5Screate_simple(hdf5_dims, count_array, 
						count_array);

	  if (mem_space_id < 0)
	    error = true;
	  else
	    {
	    if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_space_id, 
			file_space_id, H5P_DEFAULT, buffer) < 0)
	      error = true;

	    H5Sclose(mem_space_id);
	    mem_space_id = -1;
	    }
	  }

	delete[] start_array;
	start_array = 0;

	delete[] stride_array;
	stride_array = 0;
    
	delete[] count_array;
	count_array = 0;

	delete[] block_array;
	block_array = 0;

	H5Sclose(file_space_id);
	file_space_id = -1;
	}
      }

    const clock_t end_cpu_clock = clock();  // Time at end

    accum_file_read_cpu_clock += (end_cpu_clock - begin_cpu_clock);

    if (error)
      throw Except("Unable to read data", __FILE__, __LINE__);
    else
      {
      buffer_loaded = true;
      buffer_first_index = new_buffer_first_index;
      }
    }
  }


void Dataset::store_block_data(const unsigned int block_index,
			       BlockData_Double & block) const
  {
  // block must already have the correct dimensions

  if (cache_size == 0)
    {
    // Use direct write to block -- correct for Fortran ordering!

    block.set_data_reverse(block_size, get_buffer_block_ptr(block_index));
    }
  else
    {
    // Use the std::vector data retreival function

    std::vector<double> block_data(block_size);

    get_block_data(block_index, block_data);

    block.set_data(block_data);
    }
  }


const double * Dataset::get_buffer_block_ptr(const unsigned int block_index)
  const
  {
  // WARNING: The block elements are stored in Fortran-style order

  load_block(block_index);

  const unsigned int start_loc 
    = (block_index - buffer_first_index) * block_size;

  const double * array_ptr = &(buffer[start_loc]);

  return array_ptr;
  }


}  // End namspace File
}  // End namespace QuickFlash
