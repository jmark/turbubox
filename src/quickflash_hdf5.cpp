// C++ program file quickflash_hdf5.cpp

/*
  By Nathan C. Hearn
     September 1, 2006

  HDF5 utility functions.
*/


#include "quickflash_hdf5.hpp"
#include <hdf5.h>
#include <vector>
#include <map>
#include <algorithm>


namespace QuickFlash
{
namespace HDF5
{

const unsigned int Parameter_Name_Size = 80;  // READ FROM FILE INSTEAD???
const unsigned int Parameter_String_Size = 80;  // READ FROM FILE INSTEAD???

const char Parameter_Name_Label[] = "name";
const char Parameter_Value_Label[] = "value";


int get_dataset_dims(const hid_t dataset_id, std::vector<unsigned int> & dims)
  {
  int ret_val = -1;  // Error state

  if (dataset_id >= 0)
    {
    hid_t space_id = H5Dget_space(dataset_id);

    if (space_id >= 0)
      {
      int num_dims = H5Sget_simple_extent_dims(space_id, 0, 0);

      if (num_dims > 0)
	{
	hsize_t * dims_array = new hsize_t[num_dims];

	if (H5Sget_simple_extent_dims(space_id, dims_array, 0) >= 0)
	  {
	  dims.resize(num_dims);

	  for (int i = 0; i < num_dims; i++)
	    dims[i] = static_cast<unsigned int>(dims_array[i]);

	  ret_val = 0;  // OK
	  }

	delete[] dims_array;
	dims_array = 0;
	}

      H5Sclose(space_id);
      space_id = -1;
      }
    }

  return ret_val;
  }


int get_block_dataset_info(const hid_t dataset_id, unsigned int & num_blocks,
			   std::vector<unsigned int> & block_dims)
  {
  std::vector<unsigned int> dataset_dims;

  return get_block_dataset_info(dataset_id, dataset_dims, num_blocks, 
				block_dims);
  }


int get_block_dataset_info(const hid_t dataset_id, 
			   std::vector<unsigned int> & dataset_dims, 
			   unsigned int & num_blocks,
			   std::vector<unsigned int> & block_dims)
  {
  int ret_val = -1;  // Error state

  if (get_dataset_dims(dataset_id, dataset_dims) >= 0)
    {
    const unsigned int num_dataset_dims = dataset_dims.size();

    if (dataset_dims.size() > 1)
      {
      num_blocks = dataset_dims[0];  // Number of blocks

      // Determine the number of spatial dimensions -- remember Fortran order!

      block_dims.resize(0);

      for (unsigned int i = (num_dataset_dims - 1); i > 0; i--)
	{
	const unsigned int current_dim = dataset_dims[i];

	block_dims.push_back(current_dim);  // Reverse order
	}

      ret_val = 0;  // OK
      }
    }

  return ret_val;
  }


int read_attribute(const hid_t container_id, const char attribute_name[],
		   double & value)
  {
  int ret_val = -1;  // Error state

  // Locate the attribute

  hid_t attr_id = H5Aopen_name(container_id, attribute_name);

  if(attr_id >= 0)
    {
    // Get the number of elements in the data space

    hid_t space_id = H5Aget_space(attr_id);

    if (space_id >= 0)
      {
      const unsigned int num_elems 
	= static_cast<unsigned int>(H5Sget_simple_extent_npoints(space_id));

      H5Sclose(space_id);
      space_id = -1;

      if (num_elems == 1)
	{
	double read_value = 0.0;

	if (H5Aread(attr_id, H5T_NATIVE_DOUBLE, &read_value) >= 0)
	  {
	  value = read_value;
	  ret_val = 0;  // OK
	  }
	}
      }

    H5Aclose(attr_id);
    attr_id = -1;
    }

  return ret_val;
  }


int read_attribute(const hid_t container_id, const char attribute_name[],
		   int & value)
  {
  int ret_val = -1;  // Error state

  // Locate the attribute

  hid_t attr_id = H5Aopen_name(container_id, attribute_name);

  if(attr_id >= 0)
    {
    // Get the number of elements in the data space

    hid_t space_id = H5Aget_space(attr_id);

    if (space_id >= 0)
      {
      const unsigned int num_elems 
	= static_cast<unsigned int>(H5Sget_simple_extent_npoints(space_id));

      H5Sclose(space_id);
      space_id = -1;

      if (num_elems == 1)
	{
	int read_value = 0;

	if (H5Aread(attr_id, H5T_NATIVE_INT, &read_value) >= 0)
	  {
	  value = read_value;
	  ret_val = 0;  // OK
	  }
	}
      }

    H5Aclose(attr_id);
    attr_id = -1;
    }

  return ret_val;
  }


int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector<double> & data)
  {
  int ret_val = -1;  // Error state

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        const unsigned int set_size = H5Sget_simple_extent_npoints(space_id);

        H5Sclose(space_id);
        space_id = -1;

        if (set_size > 0)
          {
          double * dataset_data = new double[set_size];

          if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, dataset_data) >= 0)
            {
            data.resize(set_size);

            for (unsigned int i = 0; i < set_size; i++)
              data[i] = dataset_data[i];

            ret_val = 0;  // OK
            }

          delete[] dataset_data;
          dataset_data = 0;
          }
        }

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector<int> & data)
  {
  int ret_val = -1;  // Error state

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        const unsigned int set_size = H5Sget_simple_extent_npoints(space_id);

        H5Sclose(space_id);
        space_id = -1;

        if (set_size > 0)
          {
          int * dataset_data = new int[set_size];

          if (H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, dataset_data) >= 0)
            {
            data.resize(set_size);

            for (unsigned int i = 0; i < set_size; i++)
              data[i] = dataset_data[i];

            ret_val = 0;  // OK
            }

          delete[] dataset_data;
          dataset_data = 0;
          }
        }

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector< std::vector<double> > & data)
  {
  int ret_val = -1;  // Error state

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 2)
          {
          hsize_t dims[2];

          H5Sget_simple_extent_dims(space_id, dims, 0);

          const unsigned int num_rows = dims[0];
          const unsigned int num_cols = dims[1];

          double * dataset_data = new double[num_rows * num_cols];

          if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, dataset_data) >= 0)
            {
            // Set up the arrays

            data.resize(num_rows);

            for (unsigned int i = 0; i < num_rows; i++)
              {
              data[i].resize(num_cols);

              for (unsigned int j = 0; j < num_cols; j++)
                {
                // Assume row-major storage

                const unsigned int index = (i * num_cols) + j;

                data[i][j] = dataset_data[index];
                }
              }

            ret_val = 0;  // OK
            }

          delete[] dataset_data;
          dataset_data = 0;
          }

        H5Sclose(space_id);
        space_id = -1;
        }

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_data(const hid_t container_id, const char dataset_name[],
	      const unsigned int num_vector_elems,
	      std::vector< std::vector<double> > & data,
	      const double default_elem_value)
  {
  // NOTE: num_vector_elems corresonds to the range of the second index of data

  int ret_val = -1;  // Error state

  if ((num_vector_elems > 0) && (container_id >= 0))
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 2)
          {
          hsize_t dims[2];

          H5Sget_simple_extent_dims(space_id, dims, 0);

          const unsigned int num_rows = dims[0];
          const unsigned int num_cols = dims[1];

          double * dataset_data = new double[num_rows * num_cols];

	  const unsigned int extract_cols 
	    = std::min(num_cols, num_vector_elems);

	  const bool fill_default_values
	    = (num_vector_elems > extract_cols) ? true : false;

          if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, dataset_data) >= 0)
            {
            // Set up the arrays

            data.resize(num_rows);

            for (unsigned int i = 0; i < num_rows; i++)
              {
              data[i].resize(num_vector_elems);

              for (unsigned int j = 0; j < extract_cols; j++)
                {
                // Assume row-major storage

                const unsigned int index = (i * num_cols) + j;

                data[i][j] = dataset_data[index];
                }

	      if (fill_default_values)
		for (unsigned int j = extract_cols; j < num_vector_elems; j++)
		  data[i][j] = default_elem_value;
              }

            ret_val = 0;  // OK
            }

          delete[] dataset_data;
          dataset_data = 0;
          }

        H5Sclose(space_id);
        space_id = -1;
        }

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_data(const hid_t container_id, const char dataset_name[],
	      std::vector< std::vector<int> > & data)
  {
  int ret_val = -1;  // Error state

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 2)
          {
          hsize_t dims[2];

          H5Sget_simple_extent_dims(space_id, dims, 0);

          const unsigned int num_rows = dims[0];
          const unsigned int num_cols = dims[1];

          int * dataset_data = new int[num_rows * num_cols];

          if (H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, dataset_data) >= 0)
            {
            // Set up the arrays

            data.resize(num_rows);

            for (unsigned int i = 0; i < num_rows; i++)
              {
              data[i].resize(num_cols);

              for (unsigned int j = 0; j < num_cols; j++)
                {
                // Assume row-major storage

                const unsigned int index = (i * num_cols) + j;

                data[i][j] = dataset_data[index];
                }
              }

            ret_val = 0;  // OK
            }

          delete[] dataset_data;
          dataset_data = 0;
          }

        H5Sclose(space_id);
        space_id = -1;
        }

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_data(const hid_t container_id, const char dataset_name[],
	      const unsigned int num_vector_elems,
	      std::vector< std::vector<int> > & data,
	      const int default_elem_value)
  {
  // NOTE: num_vector_elems corresonds to the range of the second index of data

  int ret_val = -1;  // Error state

  if ((num_vector_elems > 0) && (container_id >= 0))
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 2)
          {
          hsize_t dims[2];

          H5Sget_simple_extent_dims(space_id, dims, 0);

          const unsigned int num_rows = dims[0];
          const unsigned int num_cols = dims[1];

          int * dataset_data = new int[num_rows * num_cols];

	  const unsigned int extract_cols 
	    = std::min(num_cols, num_vector_elems);

	  const bool fill_default_values 
	    = (num_vector_elems > extract_cols) ? true : false;

          if (H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, dataset_data) >= 0)
            {
            // Set up the arrays

            data.resize(num_rows);

            for (unsigned int i = 0; i < num_rows; i++)
              {
              data[i].resize(num_vector_elems);

              for (unsigned int j = 0; j < extract_cols; j++)
                {
                // Assume row-major storage

                const unsigned int index = (i * num_cols) + j;

                data[i][j] = dataset_data[index];
                }

	      if (fill_default_values)
		for (unsigned int j = extract_cols; j < num_vector_elems; j++)
		  data[i][j] = default_elem_value;
              }

            ret_val = 0;  // OK
            }

          delete[] dataset_data;
          dataset_data = 0;
          }

        H5Sclose(space_id);
        space_id = -1;
        }

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_data_stripe(const hid_t container_id, const char dataset_name[],
		     const unsigned int column_axis, 
		     const std::vector<unsigned int> & column_coords,
		     std::vector<double> & data)
  {
  int ret_val = -1;  // Error state

  const unsigned int dims = column_coords.size();

  if (column_axis < dims)
    {
    if (container_id >= 0)
      {
      // Open the database
	
      hid_t dataset_id = H5Dopen(container_id, dataset_name);

      if (dataset_id >= 0)
	{
	// Open the dataspace

	hid_t space_id = H5Dget_space(dataset_id);

	if (space_id >= 0)
	  {
	  // Determine the number of dimensions

	  const unsigned int dataset_dims 
	    = H5Sget_simple_extent_ndims(space_id);

	  if (dataset_dims == dims)
	    {
	    // Get the dimensionality of the data set

	    hsize_t * dataset_dims_array = new hsize_t[dims];

	    H5Sget_simple_extent_dims(space_id, dataset_dims_array, 0);

	    const unsigned int num_data = dataset_dims_array[column_axis];

	    bool in_bounds = true;

	    for (unsigned int i = 0; i < dims; i++)
	      if ((i != column_axis) 
		  && (column_coords[i] >= dataset_dims_array[i]))
		in_bounds = false;
	    
	    if (in_bounds)
	      {
	      // Define the region of the column to be extracted
	      // ... assume row-major storage

	      hsize_t * start_array = new hsize_t[dims];
	      hsize_t * count_array = new hsize_t[dims];

	      for (unsigned int i = 0; i < dims; i++)
		if (i != column_axis)
		  {
		  start_array[i] = static_cast<hsize_t>(column_coords[i]);
		  count_array[i] = 1;
		  }

	      start_array[column_axis] = 0;
	      count_array[column_axis] = static_cast<hsize_t>(num_data);

	      if (H5Sselect_hyperslab(space_id, H5S_SELECT_AND, start_array,
				      0, count_array, 0) >= 0)
		{
		// Create the storage space in memory

		double * data_array = new double[num_data];

		hid_t mem_space_id = H5Screate_simple(dims, count_array,
						      count_array);

		if (mem_space_id >= 0)
		  {
		  // Read and store the data

		  if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_space_id,
			      space_id, H5P_DEFAULT, data_array) >= 0)
		    {
		    data.resize(num_data);

		    for (unsigned int index = 0; index < num_data; index++)
		      data[index] = data_array[index];

		    ret_val = 0;  // OK
		    }

		  H5Sclose(mem_space_id);
		  mem_space_id = -1;
		  }

		delete[] data_array;
		data_array = 0;
		}

	      delete[] count_array;
	      count_array = 0;

	      delete[] start_array;
	      start_array = 0;
	      }

	    delete[] dataset_dims_array;
	    dataset_dims_array = 0;
	    }

	  H5Sclose(space_id);
	  space_id = -1;
	  }

	H5Dclose(dataset_id);
	dataset_id = -1;
	}
      }
    }

  return ret_val;
  }


void strip_trailing_spaces(const unsigned int str_len, char str[])
  {
  if (str_len > 0)
    {
    bool found_end = false;

    unsigned int str_ptr = str_len - 1;

    while (!found_end)
      {
      const char current_char = str[str_ptr];

      if ((current_char == '\0') || (current_char == ' '))
	str[str_ptr] = '\0';
      else
	found_end = true;

      if (str_ptr == 0)
	found_end = true;
      else
	str_ptr--;
      }
    }
  }


void strip_trailing_spaces(std::string & str)
  {
  const unsigned int str_len = str.size();

  if (str_len > 0)
    {
    bool found_end = false;

    unsigned int str_ptr = str_len - 1;

    while (!found_end)
      {
      const char current_char = str[str_ptr];

      if ((current_char != '\0') && (current_char != ' '))
	found_end = true;
      else if (str_ptr == 0)
	found_end = true;
      else
	str_ptr--;
      }

    str.resize(str_ptr + 1);
    }
  }


int read_param_table(const hid_t container_id, const char dataset_name[],
		     std::map<std::string, bool> & params)
  {
  int ret_val = -1;  // Error state

  struct Parameter_Entry
    {
    public:
      char name[Parameter_Name_Size];
      int value;
    };

  const unsigned int name_offset = HOFFSET(Parameter_Entry, name);
  const unsigned int value_offset = HOFFSET(Parameter_Entry, value);

  const unsigned int param_bytes 
    = (Parameter_Name_Size * sizeof(char)) + sizeof(int);

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 1)
          {
          hsize_t num_entries;

          H5Sget_simple_extent_dims(space_id, &num_entries, 0);

	  // Create the datatype

	  hid_t name_string_id = H5Tcopy(H5T_C_S1);
	  H5Tset_strpad(name_string_id, H5T_STR_NULLPAD);
	  H5Tset_size(name_string_id, Parameter_Name_Size);

	  // ERROR CHECKING???

	  hid_t type_id = H5Tcreate(H5T_COMPOUND, param_bytes);

	  if (type_id >= 0)
	    {
	    H5Tinsert(type_id, Parameter_Name_Label, name_offset, 
		      name_string_id);
	    H5Tinsert(type_id, Parameter_Value_Label, value_offset, 
		      H5T_NATIVE_INT);

	    Parameter_Entry * param_list = new Parameter_Entry[num_entries];

	    if (H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			param_list) >= 0)
	      {
	      for (unsigned int i = 0; i < num_entries; i++)
		{
		Parameter_Entry & entry = param_list[i];

		// Null-terminate the string and store it

		strip_trailing_spaces(Parameter_Name_Size, entry.name);

		const std::string param_name = entry.name;

		params[param_name] = (entry.value != 0) ? true : false;
		}

	      ret_val = 0;  // OK
	      }

	    delete[] param_list;
	    param_list = 0;

	    H5Tclose(type_id);
	    type_id = -1;
	    }

	  H5Tclose(name_string_id);
	  name_string_id = -1;
	  }

	H5Sclose(space_id);
	space_id = -1;
	}

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_param_table(const hid_t container_id, const char dataset_name[],
		     std::map<std::string, int> & params)
  {
  int ret_val = -1;  // Error state

  struct Parameter_Entry
    {
    public:
      char name[Parameter_Name_Size];
      int value;
    };

  const unsigned int name_offset = HOFFSET(Parameter_Entry, name);
  const unsigned int value_offset = HOFFSET(Parameter_Entry, value);

  const unsigned int param_bytes 
    = (Parameter_Name_Size * sizeof(char)) + sizeof(int);

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 1)
          {
          hsize_t num_entries;

          H5Sget_simple_extent_dims(space_id, &num_entries, 0);

	  // Create the datatype

	  hid_t name_string_id = H5Tcopy(H5T_C_S1);
	  H5Tset_strpad(name_string_id, H5T_STR_NULLPAD);
	  H5Tset_size(name_string_id, Parameter_Name_Size);

	  // ERROR CHECKING???

	  hid_t type_id = H5Tcreate(H5T_COMPOUND, param_bytes);

	  if (type_id >= 0)
	    {
	    H5Tinsert(type_id, Parameter_Name_Label, name_offset, 
		      name_string_id);
	    H5Tinsert(type_id, Parameter_Value_Label, value_offset, 
		      H5T_NATIVE_INT);

	    Parameter_Entry * param_list = new Parameter_Entry[num_entries];

	    if (H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			param_list) >= 0)
	      {
	      for (unsigned int i = 0; i < num_entries; i++)
		{
		Parameter_Entry & entry = param_list[i];

		// Null-terminate the string and store it

		strip_trailing_spaces(Parameter_Name_Size, entry.name);

		const std::string param_name = entry.name;

		params[param_name] = entry.value;
		}

	      ret_val = 0;  // OK
	      }

	    delete[] param_list;
	    param_list = 0;

	    H5Tclose(type_id);
	    type_id = -1;
	    }

	  H5Tclose(name_string_id);
	  name_string_id = -1;
	  }

	H5Sclose(space_id);
	space_id = -1;
	}

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_param_table(const hid_t container_id, const char dataset_name[],
		     std::map<std::string, double> & params)
  {
  int ret_val = -1;  // Error state

  struct Parameter_Entry
    {
    public:
      char name[Parameter_Name_Size];
      double value;
    };

  const unsigned int name_offset = HOFFSET(Parameter_Entry, name);
  const unsigned int value_offset = HOFFSET(Parameter_Entry, value);

  const unsigned int param_bytes 
    = (Parameter_Name_Size * sizeof(char)) + sizeof(double);

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 1)
          {
          hsize_t num_entries;

          H5Sget_simple_extent_dims(space_id, &num_entries, 0);

	  // Create the datatype

	  hid_t name_string_id = H5Tcopy(H5T_C_S1);
	  H5Tset_strpad(name_string_id, H5T_STR_NULLPAD);
	  H5Tset_size(name_string_id, Parameter_Name_Size);

	  // ERROR CHECKING???

	  hid_t type_id = H5Tcreate(H5T_COMPOUND, param_bytes);

	  if (type_id >= 0)
	    {
	    H5Tinsert(type_id, Parameter_Name_Label, name_offset, 
		      name_string_id);
	    H5Tinsert(type_id, Parameter_Value_Label, value_offset, 
		      H5T_NATIVE_DOUBLE);

	    Parameter_Entry * param_list = new Parameter_Entry[num_entries];

	    if (H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			param_list) >= 0)
	      {
	      for (unsigned int i = 0; i < num_entries; i++)
		{
		Parameter_Entry & entry = param_list[i];

		// Null-terminate the string and store it

		strip_trailing_spaces(Parameter_Name_Size, entry.name);

		const std::string param_name = entry.name;

		params[param_name] = entry.value;
		}

	      ret_val = 0;  // OK
	      }

	    delete[] param_list;
	    param_list = 0;

	    H5Tclose(type_id);
	    type_id = -1;
	    }

	  H5Tclose(name_string_id);
	  name_string_id = -1;
	  }

	H5Sclose(space_id);
	space_id = -1;
	}

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_param_table(const hid_t container_id, const char dataset_name[],
		     std::map<std::string, std::string> & params)
  {
  int ret_val = -1;  // Error state

  struct Parameter_Entry
    {
    public:
      char name[Parameter_Name_Size];
      char value[Parameter_String_Size];
    };

  const unsigned int name_offset = HOFFSET(Parameter_Entry, name);
  const unsigned int value_offset = HOFFSET(Parameter_Entry, value);

  const unsigned int param_bytes 
    = (Parameter_Name_Size * sizeof(char)) 
      + (Parameter_String_Size * sizeof(char));

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 1)
          {
          hsize_t num_entries;

          H5Sget_simple_extent_dims(space_id, &num_entries, 0);

	  // Create the datatype

	  hid_t name_string_id = H5Tcopy(H5T_C_S1);
	  H5Tset_strpad(name_string_id, H5T_STR_NULLPAD);
	  H5Tset_size(name_string_id, Parameter_Name_Size);

	  hid_t value_string_id = H5Tcopy(H5T_C_S1);
	  H5Tset_strpad(value_string_id, H5T_STR_NULLPAD);
	  H5Tset_size(value_string_id, Parameter_String_Size);

	  // ERROR CHECKING???

	  hid_t type_id = H5Tcreate(H5T_COMPOUND, param_bytes);

	  if (type_id >= 0)
	    {
	    H5Tinsert(type_id, Parameter_Name_Label, name_offset, 
		      name_string_id);
	    H5Tinsert(type_id, Parameter_Value_Label, value_offset, 
		      value_string_id);

	    Parameter_Entry * param_list = new Parameter_Entry[num_entries];

	    if (H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			param_list) >= 0)
	      {
	      for (unsigned int i = 0; i < num_entries; i++)
		{
		Parameter_Entry & entry = param_list[i];

		// Null-terminate the string and store it

		strip_trailing_spaces(Parameter_Name_Size, entry.name);
		strip_trailing_spaces(Parameter_String_Size, entry.value);

		const std::string param_name = entry.name;

		params[param_name] = entry.value;
		}

	      ret_val = 0;  // OK
	      }

	    delete[] param_list;
	    param_list = 0;

	    H5Tclose(type_id);
	    type_id = -1;
	    }

	  H5Tclose(name_string_id);
	  name_string_id = -1;
	  }

	H5Sclose(space_id);
	space_id = -1;
	}

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_compound_param(const hid_t container_id, const char dataset_name[],
			const char elem_name[], bool & value)
  {
  int ret_val = -1;  // Error state

  struct Parameter_Entry
    {
    public:
      int value;
    };

  const unsigned int value_offset = HOFFSET(Parameter_Entry, value);

  const unsigned int param_bytes = sizeof(int);

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 1)
          {
          hsize_t num_entries;

          H5Sget_simple_extent_dims(space_id, &num_entries, 0);

	  if (num_entries == 1)
	    {
	    // Create the datatype

	    hid_t type_id = H5Tcreate(H5T_COMPOUND, param_bytes);

	    if (type_id >= 0)
	      {
	      H5Tinsert(type_id, elem_name, value_offset, H5T_NATIVE_INT);

	      Parameter_Entry entry;

	      if (H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			  &entry) >= 0)
		{
		value = (entry.value != 0) ? true : false;

		ret_val = 0;  // OK
		}

	      H5Tclose(type_id);
	      type_id = -1;
	      }
	    }
	  }

	H5Sclose(space_id);
	space_id = -1;
	}

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_compound_param(const hid_t container_id, const char dataset_name[],
			const char elem_name[], int & value)
  {
  int ret_val = -1;  // Error state

  struct Parameter_Entry
    {
    public:
      int value;
    };

  const unsigned int value_offset = HOFFSET(Parameter_Entry, value);

  const unsigned int param_bytes = sizeof(int);

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 1)
          {
          hsize_t num_entries;

          H5Sget_simple_extent_dims(space_id, &num_entries, 0);

	  if (num_entries == 1)
	    {
	    // Create the datatype

	    hid_t type_id = H5Tcreate(H5T_COMPOUND, param_bytes);

	    if (type_id >= 0)
	      {
	      H5Tinsert(type_id, elem_name, value_offset, H5T_NATIVE_INT);

	      Parameter_Entry entry;

	      if (H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			  &entry) >= 0)
		{
		value = entry.value;

		ret_val = 0;  // OK
		}

	      H5Tclose(type_id);
	      type_id = -1;
	      }
	    }
	  }

	H5Sclose(space_id);
	space_id = -1;
	}

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_compound_param(const hid_t container_id, const char dataset_name[],
			const char elem_name[], double & value)
  {
  int ret_val = -1;  // Error state

  struct Parameter_Entry
    {
    public:
      double value;
    };

  const unsigned int value_offset = HOFFSET(Parameter_Entry, value);

  const unsigned int param_bytes = sizeof(double);

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 1)
          {
          hsize_t num_entries;

          H5Sget_simple_extent_dims(space_id, &num_entries, 0);

	  if (num_entries == 1)
	    {
	    // Create the datatype

	    hid_t type_id = H5Tcreate(H5T_COMPOUND, param_bytes);

	    if (type_id >= 0)
	      {
	      H5Tinsert(type_id, elem_name, value_offset, H5T_NATIVE_DOUBLE);

	      Parameter_Entry entry;

	      if (H5Dread(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			  &entry) >= 0)
		{
		value = entry.value;

		ret_val = 0;  // OK
		}

	      H5Tclose(type_id);
	      type_id = -1;
	      }
	    }
	  }

	H5Sclose(space_id);
	space_id = -1;
	}

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int read_string_set(const hid_t container_id, const char dataset_name[],
		    std::vector<std::string> & data)
  {
  int ret_val = -1;  // Error state

  if (container_id >= 0)
    {
    // Open the database

    hid_t dataset_id = H5Dopen(container_id, dataset_name);

    if (dataset_id >= 0)
      {
      // Open the dataspace

      hid_t space_id = H5Dget_space(dataset_id);

      if (space_id >= 0)
        {
        // Determine the number of dimensions

        const int num_dims = H5Sget_simple_extent_ndims(space_id);

        if (num_dims == 2)
          {
          hsize_t dims[2];

          H5Sget_simple_extent_dims(space_id, dims, 0);

	  const unsigned int num_entries = dims[0];
	  const unsigned int entry_elems = dims[1];

	  if (entry_elems == 1)
	    {
	    // Find the maximum string length

	    hid_t type_id = H5Dget_type(dataset_id);

	    if (type_id >= 0)
	      {
              // Add an extra character for null (needed for some compilers)

	      const unsigned int max_str_len = H5Tget_size(type_id) + 1;

	      hid_t local_type_id = H5Tcopy(H5T_C_S1);

	      if (local_type_id >= 0)
		{
		H5Tset_size(local_type_id, max_str_len);

		char * dataset_data = new char[num_entries * max_str_len];

		data.resize(num_entries);

		if (H5Dread(dataset_id, local_type_id, H5S_ALL, H5S_ALL,
			    H5P_DEFAULT, dataset_data) >= 0)
		  {
		  for (unsigned int i = 0; i < num_entries; i++)
		    {
		    const unsigned int start_index = i * max_str_len;

		    const char * char_ptr = &(dataset_data[start_index]);

		    data[i].assign(char_ptr, max_str_len);

		    strip_trailing_spaces(data[i]);
		    }

		  ret_val = 0;  // OK
		  }

		delete[] dataset_data;
		dataset_data = 0;

		H5Tclose(local_type_id);
		local_type_id = -1;
		}

	      H5Tclose(type_id);
	      type_id = -1;
	      }
	    }
	  }

        H5Sclose(space_id);
        space_id = -1;
        }

      H5Dclose(dataset_id);
      dataset_id = -1;
      }
    }

  return ret_val;
  }


int get_num_members(const hid_t container_id)
  {
  int ret_val = -1;  // Error state

  if (container_id >= 0)
    {
    hsize_t num_members = 0;

    if (H5Gget_num_objs(container_id, &num_members) >= 0)
      {
      ret_val = static_cast<int>(num_members);
      }
    }

  return ret_val;
  }


int member_index(const hid_t container_id, const std::string & name)
  {
  int ret_val = -1;  // Error state

  if (container_id >= 0)
    {
    const int num_members = get_num_members(container_id);

    if (num_members > 0)
      {
      int member_index = 0;

      bool present = false;

      while ((member_index < num_members) && (!present))
	{
	const int string_length 
	  = H5Gget_objname_by_idx(container_id, member_index, 0, 0) + 1;

	if (string_length > 0)
	  {
	  char * string_container = new char[string_length];

	  H5Gget_objname_by_idx(container_id, member_index, string_container,
				string_length);

	  if (name == string_container)
	    {
	    present = true;
	    ret_val = member_index;
	    }

	  delete[] string_container;
	  string_container = 0;
	  }

	member_index++;
	}
      }
    }

  return ret_val;
  }


int member_index(const hid_t container_id, const char name[])
  {
  const std::string name_str = name;

  return member_index(container_id, name_str);
  }


bool member_present(const hid_t container_id, const std::string & name)
  {
  bool ret_val = false;  // Default state

  if (member_index(container_id, name) >= 0)
    ret_val = true;

  return ret_val;
  }


bool member_present(const hid_t container_id, const char name[])
  {
  const std::string name_str = name;

  return member_present(container_id, name_str);
  }


}  // End namespace HDF5
}  // End namespace QuickFlash
