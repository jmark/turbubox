// C++ program file quickflash_stats_bins.cpp

/*
  By Nathan C. Hearn
     February 18, 2008

  Statistics routines for binned data.
*/


#include "quickflash_stats_bins.hpp"
#include <cmath>
#include <vector>
#include "quickflash_except.hpp"
#include "quickflash_geometry.hpp"
#include "quickflash_block_utils.hpp"
#include "quickflash_counters.hpp"


namespace QuickFlash
{
namespace Stats
{

// Class Bin

void Bin::reset()
  {
  num_entries = 0;

  sum_weight = 0.0;

  sum_weight_values = 0.0;
  sum_weight_sqValues = 0.0;

  min_value = 0.0;
  max_value = 0.0;
  }


void Bin::reset(const Bin & source)
  {
  if (&source != this)
    {
    num_entries = source.num_entries;

    sum_weight = source.sum_weight;

    sum_weight_values = source.sum_weight_values;
    sum_weight_sqValues = source.sum_weight_sqValues;
    }
  }


void Bin::add_bin(const Bin & source)
  {
  num_entries += source.num_entries;

  sum_weight += source.sum_weight;

  sum_weight_values += source.sum_weight_values;
  sum_weight_sqValues += source.sum_weight_sqValues;

  const double source_min_val = source.min_value;
  const double source_max_val = source.max_value;

  if (source_min_val < min_value)
    min_value = source_min_val;

  if (source_max_val > max_value)
    max_value = source_max_val;
  }


void Bin::add_entry(const double value, const double weight)
  {
  // SHOULD WE ALLOW ZERO (OR NEGATIVE) WEIGHTS???

  sum_weight += weight;

  const double weight_value = weight * value;

  sum_weight_values += weight_value;

  sum_weight_sqValues += weight_value * value;

  if (num_entries > 0)
    {
    if (value < min_value)
      min_value = value;
    else if (value > max_value)
      max_value = value;
    }
  else
    {
    min_value = value;
    max_value = value;
    }

  num_entries++;
  }


double Bin::get_pop_variance() const
  {
  double variance = 0.0;

  if (sum_weight > 0.0)
    {
    const double oneOver_sum_weight = 1.0 / sum_weight;

    const double mean_val = get_mean();

    variance = (sum_weight_sqValues * oneOver_sum_weight) 
               - (mean_val * mean_val);

    if (variance < 0.0)
      variance = 0.0;
    }

  return variance;
  }


double Bin::get_pop_mean_square_difference(const double reference_value) const
  {
  // With accumulated quantities, the mean square difference from reference
  // value q is given by
  //
  //      (1 / W) sum_i[w_i (x_i - q)^2]
  //          = [(1 / W) sum_i(w_i x_i^2)] - 2 q <x> + q^2
  //
  // where <x> is the weighted mean value of x_i

  double mean_square = 0.0;

  if (sum_weight > 0.0)
    {
    const double oneOver_sum_weight = 1.0 / sum_weight;

    const double mean_val = get_mean();

    mean_square = (sum_weight_sqValues * oneOver_sum_weight) 
                  + (reference_value * (reference_value - (2.0 * mean_val)));

    if (mean_square < 0.0)
      mean_square = 0.0;
    }

  return mean_square;
  }


// Class BinArray

void BinArray::reset()
  {
  bin_count = 0;

  min_array_bounds = 0.0;
  max_array_bounds = 0.0;

  bin_min_bounds.clear();

  bin_center_coords.clear();
  bin_widths.clear();

  bins.clear();
  }


void BinArray::reset(const double min_boundary, const unsigned int num_bins,
		     const double bin_width)
  {
  if (num_bins < 1)
    throw Except("No bins requested", __FILE__, __LINE__);

  if (!(bin_width > 0.0))
    throw Except("Bin widths must be positive", __FILE__, __LINE__);

  // Build the arrays of boundaries

  std::vector<double> bin_min_boundary_list(num_bins);

  for (unsigned int index = 0; index < num_bins; index++)
    bin_min_boundary_list[index] 
      = min_boundary + (static_cast<double>(index) * bin_width);

  const double max_boundary 
    = min_boundary + (static_cast<double>(num_bins) * bin_width);

  reset(bin_min_boundary_list, max_boundary);
  }


void BinArray::reset(const double min_boundary, const double max_boundary, 
		     const unsigned int num_bins, const bool log_scale)
  {
  if (num_bins < 1)
    throw Except("No bins requested", __FILE__, __LINE__);

  if (!(min_boundary < max_boundary))
    throw Except("Min boundary must be less than max boundary", __FILE__,
		 __LINE__);

  if (log_scale && (!(min_boundary > 0.0)))
    throw Except("Array boundaries must be positive for log scale", __FILE__, 
		 __LINE__);

  // Find (the logarithms of) the start and stop positions, and the bin widths

  double array_start = 0.0;
  double array_end = 0.0;
  double array_width = 0.0;

  if (log_scale)
    {
    array_start = std::log(min_boundary);
    array_end = std::log(max_boundary);
    }
  else
    {
    array_start = min_boundary;
    array_end = max_boundary;
    }

  array_width = array_end - array_start;

  if (!(array_width > 0.0))
    throw Except("Unable find bin size", __FILE__, __LINE__);

  const double bin_width = array_width / static_cast<double>(num_bins);

  // Get the minimum bin boundary coords

  std::vector<double> bin_min_boundary_list(num_bins);

  for (unsigned int index = 0; index < num_bins; index++)
    {
    double current_min_boundary 
      = array_start + (static_cast<double>(index) * bin_width);

    if (log_scale)
      current_min_boundary = std::exp(current_min_boundary);

    bin_min_boundary_list[index] = current_min_boundary;
    }

  reset(bin_min_boundary_list, max_boundary);
  }


void BinArray::reset(const std::vector<double> bin_min_boundary_list, 
		     const double max_boundary)
  {
  const unsigned int num_bins = bin_min_boundary_list.size();

  if (num_bins < 1)
    throw Except("No bins specified", __FILE__, __LINE__);

  // Check the boundaries

  double current_bin_min = bin_min_boundary_list[0];

  for (unsigned int index = 1; index < num_bins; index++)
    {
    const double next_bin_min = bin_min_boundary_list[index];

    if (!(next_bin_min > current_bin_min))
      throw Except("Incompatible bin min boundaries found", __FILE__, 
		   __LINE__);

    current_bin_min = next_bin_min;
    }

  if (!(max_boundary > bin_min_boundary_list[num_bins - 1]))
    throw Except("Max boundary not greater than last bin min", __FILE__, 
		 __LINE__);

  // Set up the bins
  
  bin_count = num_bins;

  min_array_bounds = bin_min_boundary_list[0];
  max_array_bounds = max_boundary;

  bin_min_bounds.resize(num_bins);

  bin_center_coords.resize(num_bins);
  bin_widths.resize(num_bins);

  bins.resize(num_bins);

  const unsigned int last_bin_index = num_bins - 1;

  for (unsigned int index = 0; index < last_bin_index; index++)
    {
    const double bin_min_coord = bin_min_boundary_list[index];
    const double bin_max_coord = bin_min_boundary_list[index + 1];

    bin_min_bounds[index] = bin_min_coord;

    bin_center_coords[index] = 0.5 * (bin_min_coord + bin_max_coord);
    bin_widths[index] = bin_max_coord - bin_min_coord;

    bins[index].reset();
    }

  const double last_bin_min = bin_min_boundary_list[last_bin_index];

  bin_min_bounds[last_bin_index] = last_bin_min;

  bin_center_coords[last_bin_index] = 0.5 * (last_bin_min + max_boundary);
  bin_widths[last_bin_index] = max_boundary - last_bin_min;

  bins[last_bin_index].reset();
  }


void BinArray::reset(const BinArray & source)
  {
  if (&source != this)
    {
    bin_count = source.bin_count;

    min_array_bounds = source.min_array_bounds;
    max_array_bounds = source.max_array_bounds;

    bin_min_bounds = source.bin_min_bounds;

    bin_center_coords = source.bin_center_coords;
    bin_widths = source.bin_widths;

    bins = source.bins;
    }
  }


void BinArray::reset_bins()
  {
  // Reset bin values but leave array intact

  for (unsigned int bin_index = 0; bin_index < bin_count; bin_index++)
    bins[bin_index].reset();
  }


void BinArray::add_data(const double coordinate, const double value,
			const double weight)
  {
  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  if (in_bounds(coordinate))
    {
    const unsigned int bin_index = get_bin_index(coordinate);

    bins[bin_index].add_entry(value, weight);
    }
  }


void BinArray::add_data_spread(const double coordinate, const double value,
			       const double spread, const double weight)
  {
  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  if (!(spread > 0.0))
    add_data(coordinate, value, weight);
  else
    {
    const double half_spread = 0.5 * spread;

    const double min_coord = coordinate - half_spread;
    const double max_coord = coordinate + half_spread;

    if ((!(max_coord < min_array_bounds)) && (min_coord < max_array_bounds))
      {
      const unsigned int min_bin_index = get_nearest_bin_index(min_coord);
      const unsigned int max_bin_index = get_nearest_bin_index(max_coord);

      const double oneOver_spread = 1.0 / spread;

      for (unsigned int index = min_bin_index; index <= max_bin_index; index++)
	{
	// See how much of the spread is in each bin and modify the weight

	const double overlap = get_overlap_area(index, coordinate, spread);

	const double bin_weight = (overlap * oneOver_spread) * weight;

	bins[index].add_entry(value, bin_weight);
	}
      }
    }
  }


bool BinArray::in_bounds(const double coordinate, const double spread) 
  const
  {
  bool ret_val = false;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  if (!(spread > 0.0))
    {
    if ((!(coordinate < min_array_bounds)) && (coordinate < max_array_bounds))
      ret_val = true;
    }
  else
    {
    const double half_spread = 0.5 * spread;

    const double min_coord = coordinate - half_spread;
    const double max_coord = coordinate + half_spread;

    if ((!(max_coord < min_array_bounds)) && (min_coord < max_array_bounds))
      ret_val = true;
    }

  return ret_val;
  }


unsigned int BinArray::get_bin_index(const double coordinate) const
  {
  unsigned int ret_val = bin_count;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  if (!in_bounds(coordinate))
    throw Except("Coordinate out of bounds", __FILE__, __LINE__);

  // Find the bin containing coordinate

  const unsigned int last_bin_index = bin_count - 1;

  unsigned int lo_index = 0;
  unsigned int hi_index = last_bin_index;

  double lo_bound = min_array_bounds;
  double hi_bound = max_array_bounds;

  while (lo_index < (hi_index - 1))
    {
    const unsigned int center_index = (lo_index + hi_index) / 2;

    const double bin_min = bin_min_bounds[center_index];

    if (coordinate < bin_min)
      {
      hi_index = center_index;
      hi_bound = bin_min;
      }
    else
      {
      lo_index = center_index;
      lo_bound = bin_min;
      }
    }

  // If two bins left in range, determine which contains coordinate

  if (lo_index != hi_index)
    {
    if (!(coordinate < bin_min_bounds[hi_index]))
      ret_val = hi_index;
    else
      ret_val = lo_index;
    }
  else
    ret_val = lo_index;

  return ret_val;
  }


unsigned int BinArray::get_nearest_bin_index(const double coordinate) const
  {
  unsigned int ret_val = 0;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  if (coordinate < min_array_bounds)
    ret_val = 0;
  else if (!(coordinate < max_array_bounds))
    ret_val = bin_count - 1;
  else
    ret_val = get_bin_index(coordinate);

  return ret_val;
  }


bool BinArray::in_bin(const unsigned int bin_index, const double coordinate, 
		      const double spread) const
  {
  bool ret_val = false;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  const double bin_min = bin_min_bounds[bin_index];

  double bin_max = 0.0;

  if (bin_index == (bin_count - 1))
    bin_max = max_array_bounds;
  else
    bin_max = bin_min_bounds[bin_index + 1];

  if (!(spread > 0.0))
    {
    if ((!(coordinate < bin_min)) && (coordinate < bin_max))
      ret_val = true;
    }
  else
    {
    const double half_spread = 0.5 * spread;

    const double min_coord = coordinate - half_spread;
    const double max_coord = coordinate + half_spread;

    if ((!(max_coord < bin_min)) && (min_coord < bin_max))
      ret_val = true;
    }

  return ret_val;
  }


double BinArray::get_overlap_area(const unsigned int bin_index, 
				  const double coordinate,
				  const double spread) const
  {
  double ret_val = 0.0;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  if (spread > 0.0)
    {
    const double bin_min = bin_min_bounds[bin_index];

    double bin_max = 0.0;

    if (bin_index == (bin_count - 1))
      bin_max = max_array_bounds;
    else
      bin_max = bin_min_bounds[bin_index + 1];

    const double half_spread = 0.5 * spread;

    double min_overlap_coord = coordinate - half_spread;
    double max_overlap_coord = coordinate + half_spread;

    if (min_overlap_coord < bin_min)
      min_overlap_coord = bin_min;

    if (max_overlap_coord > bin_max)
      max_overlap_coord = bin_max;

    const double overlap_area = max_overlap_coord - min_overlap_coord;

    if (overlap_area > 0.0)
      ret_val = overlap_area;
    }

  return ret_val;
  }


double BinArray::get_peak_mean(unsigned int & bin_index) const
  {
  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  // Find the bin with the largest mean value (ignoring bins without data)

  bool data_found = false;

  double value_max = 0.0;
  unsigned int index_max = 0;

  for (unsigned int index = 0; index < bin_count; index++)
    {
    const Bin & current_bin = bins[index];

    if (current_bin.get_weight_sum() > 0.0)
      {
      const double bin_mean = current_bin.get_mean();

      if (data_found)
	{
	if (bin_mean > value_max)
	  {
	  value_max = bin_mean;
	  index_max = index;
	  }
	}
      else
	{
	value_max = bin_mean;
	index_max = index;
	
	data_found = true;
	}
      }
    }

  if (!data_found)
    throw Except("No data found in array", __FILE__, __LINE__);

  bin_index = index_max;

  return value_max;
  }


double BinArray::get_transition_pos(const double threshold_value,
				    const bool start_lower_bound) const
  {
  double ret_val = 0.0;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  // Find the first bin with data (positive weight)

  unsigned int start_index = 0;

  bool found_start = false;

  if (start_lower_bound)
    {
    while ((start_index < bin_count) && (!found_start))
      {
      if (bins[start_index].get_weight_sum() > 0.0)
	found_start = true;
      else
	start_index++;
      }
    }
  else
    {
    start_index = bin_count;  // Points to bin after current bin

    while ((start_index > 0) && (!found_start))
      {
      if (bins[start_index - 1].get_weight_sum() > 0.0)
	found_start = true;
      else
	start_index --;
      }

    start_index--;  // Move to current bin
    }

  if (!found_start)
    throw Except("Unable to find starting bin for comparison", __FILE__, 
		 __LINE__);

  const double start_value = bins[start_index].get_mean();

  const double init_delta_value = threshold_value - start_value;

  if (!(std::fabs(init_delta_value) > 0.0))
    ret_val = bin_center_coords[start_index];  // At transition point ...
  else
    {
    unsigned int lo_index = 0;
    unsigned int hi_index = 0;

    bool found_transition = false;

    if (start_lower_bound)
      {
      unsigned int bin_index = start_index + 1;
      unsigned int last_valid_bin = start_index;

      while ((bin_index < bin_count) && (!found_transition))
	{
	const Bin & current_bin = bins[bin_index];

	if (current_bin.get_weight_sum() > 0.0)
	  {
	  const double current_value = current_bin.get_mean();

	  const double current_delta_value = threshold_value - current_value;

	  if (!((current_delta_value * init_delta_value) > 0.0))
	    found_transition = true;
	  else
	    {
	    last_valid_bin = bin_index;
	    bin_index++;
	    }
	  }
	else
	  bin_index++;
	}

      if (found_transition)
	{
	lo_index = last_valid_bin;
	hi_index = bin_index;
	}
      }
    else
      {
      unsigned int bin_index = start_index;  // One after current_bin
      unsigned int last_valid_bin = start_index;

      while ((bin_index > 0) && (!found_transition))
	{
	bin_index--;

	const Bin & current_bin = bins[bin_index];

	if (current_bin.get_weight_sum() > 0.0)
	  {
	  const double current_value = current_bin.get_mean();

	  const double current_delta_value = threshold_value - current_value;

	  if (!((current_delta_value * init_delta_value) > 0.0))
	    found_transition = true;
	  else
	    last_valid_bin = bin_index;
	  }
	}

      if (found_transition)
	{
	lo_index = bin_index;
	hi_index = last_valid_bin;
	}
      }

    if (!found_transition)
      throw Except("Unable to find transition", __FILE__, __LINE__);
    else
      {
      const double lo_pos = bin_center_coords[lo_index];
      const double hi_pos = bin_center_coords[hi_index];

      const double lo_val = bins[lo_index].get_mean();
      const double hi_val = bins[hi_index].get_mean();

      ret_val = find_threshold_pos(lo_pos, hi_pos, lo_val, hi_val, 
				   threshold_value);
      }
    }
    
  return ret_val;
  }


double BinArray::get_first_pos_below(const double threshold_value,
				     const bool start_lower_bound) const
  {
  double ret_val = 0.0;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  // Find the first bin with data (positive weight)

  unsigned int start_index = 0;

  bool found_start = false;

  if (start_lower_bound)
    {
    while ((start_index < bin_count) && (!found_start))
      {
      if (bins[start_index].get_weight_sum() > 0.0)
	found_start = true;
      else
	start_index++;
      }
    }
  else
    {
    start_index = bin_count;  // Points to bin after current bin

    while ((start_index > 0) && (!found_start))
      {
      if (bins[start_index - 1].get_weight_sum() > 0.0)
	found_start = true;
      else
	start_index --;
      }

    start_index--;  // Move to current bin
    }

  if (!found_start)
    throw Except("Unable to find starting bin for comparison", __FILE__, 
		 __LINE__);

  const double start_value = bins[start_index].get_mean();

  if (!(start_value > threshold_value))
    ret_val = bin_center_coords[start_index];  // At transition point ...
  else
    {
    unsigned int lo_index = 0;
    unsigned int hi_index = 0;

    bool found_first_pos = false;

    if (start_lower_bound)
      {
      unsigned int bin_index = start_index + 1;
      unsigned int last_valid_bin = start_index;

      while((bin_index < bin_count) && (!found_first_pos))
	{
	const Bin & current_bin = bins[bin_index];

	if (current_bin.get_weight_sum() > 0.0)
	  {
	  const double current_value = current_bin.get_mean();

	  if (!(current_value > threshold_value))
	    found_first_pos = true;
	  else
	    {
	      last_valid_bin = bin_index;
	      bin_index++;
	    }
	  }
	else
	  bin_index++;
	}

      if (found_first_pos)
	{
	lo_index = last_valid_bin;
	hi_index = bin_index;
	}
      }
    else
      {
      unsigned int bin_index = start_index;  // One after current_bin
      unsigned int last_valid_bin = start_index;

      while ((bin_index > 0) && (!found_first_pos))
	{
	bin_index--;

	const Bin & current_bin = bins[bin_index];

	if (current_bin.get_weight_sum() > 0.0)
	  {
	  const double current_value = current_bin.get_mean();

	  if (!(current_value > threshold_value))
	    found_first_pos = true;
	  else
	    last_valid_bin = bin_index;
	  }
	}

      if (found_first_pos)
	{
	lo_index = bin_index;
	hi_index = last_valid_bin;
	}
      }

    if (!found_first_pos)
      throw Except("Unable to find transition", __FILE__, __LINE__);
    else
      {
      const double lo_pos = bin_center_coords[lo_index];
      const double hi_pos = bin_center_coords[hi_index];

      const double lo_val = bins[lo_index].get_mean();
      const double hi_val = bins[hi_index].get_mean();

      ret_val = find_threshold_pos(lo_pos, hi_pos, lo_val, hi_val, 
				   threshold_value);
      }
    }
    
  return ret_val;
  }


double BinArray::get_first_pos_above(const double threshold_value,
				     const bool start_lower_bound) const
  {
  double ret_val = 0.0;

  if (bin_count < 1)
    throw Except("Array not set up", __FILE__, __LINE__);

  // Find the first bin with data (positive weight)

  unsigned int start_index = 0;

  bool found_start = false;

  if (start_lower_bound)
    {
    while ((start_index < bin_count) && (!found_start))
      {
      if (bins[start_index].get_weight_sum() > 0.0)
	found_start = true;
      else
	start_index++;
      }
    }
  else
    {
    start_index = bin_count;  // Points to bin after current bin

    while ((start_index > 0) && (!found_start))
      {
      if (bins[start_index - 1].get_weight_sum() > 0.0)
	found_start = true;
      else
	start_index --;
      }

    start_index--;  // Move to current bin
    }

  if (!found_start)
    throw Except("Unable to find starting bin for comparison", __FILE__, 
		 __LINE__);

  const double start_value = bins[start_index].get_mean();

  if (!(start_value < threshold_value))
    ret_val = bin_center_coords[start_index];  // At transition point ...
  else
    {
    unsigned int lo_index = 0;
    unsigned int hi_index = 0;

    bool found_first_pos = false;

    if (start_lower_bound)
      {
      unsigned int bin_index = start_index + 1;
      unsigned int last_valid_bin = start_index;

      while((bin_index < bin_count) && (!found_first_pos))
	{
	const Bin & current_bin = bins[bin_index];

	if (current_bin.get_weight_sum() > 0.0)
	  {
	  const double current_value = current_bin.get_mean();

	  if (!(current_value < threshold_value))
	    found_first_pos = true;
	  else
	    {
	      last_valid_bin = bin_index;
	      bin_index++;
	    }
	  }
	else
	  bin_index++;
	}

      if (found_first_pos)
	{
	lo_index = last_valid_bin;
	hi_index = bin_index;
	}
      }
    else
      {
      unsigned int bin_index = start_index;  // One after current_bin
      unsigned int last_valid_bin = start_index;

      while ((bin_index > 0) && (!found_first_pos))
	{
	bin_index--;

	const Bin & current_bin = bins[bin_index];

	if (current_bin.get_weight_sum() > 0.0)
	  {
	  const double current_value = current_bin.get_mean();

	  if (!(current_value < threshold_value))
	    found_first_pos = true;
	  else
	    last_valid_bin = bin_index;
	  }
	}

      if (found_first_pos)
	{
	lo_index = bin_index;
	hi_index = last_valid_bin;
	}
      }

    if (!found_first_pos)
      throw Except("Unable to find transition", __FILE__, __LINE__);
    else
      {
      const double lo_pos = bin_center_coords[lo_index];
      const double hi_pos = bin_center_coords[hi_index];

      const double lo_val = bins[lo_index].get_mean();
      const double hi_val = bins[hi_index].get_mean();

      ret_val = find_threshold_pos(lo_pos, hi_pos, lo_val, hi_val, 
				   threshold_value);
      }
    }
    
  return ret_val;
  }


double BinArray::find_threshold_pos(const double lo_pos, const double hi_pos,
				    const double lo_val, const double hi_val,
				    const double threshold_value)
  {
  double ret_val = 0.0;

  if (!(lo_pos < hi_pos))
    throw Except("Low position must be less than high position", __FILE__,
		 __LINE__);

  if (((lo_val > threshold_value) && (hi_val > threshold_value))
      || ((lo_val < threshold_value) && (hi_val < threshold_value)))
    throw Except("Boundary values do not cross threshold", __FILE__, __LINE__);


  const double delta_threshold_lo = threshold_value - lo_val;

  const double delta_pos_lo
    = delta_threshold_lo * ((hi_pos - lo_pos) / (hi_val - lo_val));

  ret_val = lo_pos + delta_pos_lo;

  return ret_val;
  }


// Class BinGrid

BinGrid::BinGrid() : dims(0), grid_info(), grid_data() { }


BinGrid::BinGrid(const std::vector<double> & grid_min_coords, 
		 const std::vector<double> & bin_grid_spacing, 
		 const std::vector<unsigned int> & bin_grid_dims) :
  dims(0), grid_info(), grid_data()
  { reset(grid_min_coords, bin_grid_spacing, bin_grid_dims); }


BinGrid::BinGrid(const BinGrid & source) :
  dims(0), grid_info(), grid_data()
  { reset(source); }


void BinGrid::reset()
  {
  dims = 0;

  grid_info.reset();
  grid_data.reset();
  }


void BinGrid::reset(const std::vector<double> & grid_min_coords, 
		    const std::vector<double> & bin_grid_spacing, 
		    const std::vector<unsigned int> & bin_grid_dims)
  {
  reset();

  dims = grid_min_coords.size();

  if (dims > 0)
    {
    if (bin_grid_spacing.size() != dims)
      throw Except("Incompatible grid spacing vector size", __FILE__, 
		   __LINE__);

    for (unsigned int i = 0; i < dims; i++)
      if (!(bin_grid_spacing[i] > 0.0))
	throw Except("Grid spacing elements must be positive", __FILE__,
		     __LINE__);

    if (bin_grid_dims.size() != dims)
      throw Except("Incompatible grid dimensions vector size", __FILE__, 
		   __LINE__);

    for (unsigned int i = 0; i < dims; i++)
      if (bin_grid_dims[i] < 1)
	throw Except("Grid dimensions elements must be positive", __FILE__,
		     __LINE__);

    std::vector<double> grid_max_coords = grid_min_coords;

    for (unsigned int i = 0; i < dims; i++)
      {
      const double axis_grid_width 
	= bin_grid_spacing[i] * static_cast<double>(bin_grid_dims[i]);

      grid_max_coords[i] += axis_grid_width;
      }

    grid_info.reset(grid_min_coords, grid_max_coords, bin_grid_dims,
		    Geometry::Cartesian);

    grid_data.reset(bin_grid_dims);
    }
  }


void BinGrid::reset(const BinGrid & source)
  {
  if (&source != this)
    {
    dims = source.dims;

    grid_info = source.grid_info;
    grid_data = source.grid_data;
    }
  }


void BinGrid::add_data_spread(const std::vector<double> & data_center_pos,
			      const double data_value,
			      const std::vector<double> & data_spread,
			      const double weight)
  {
  if (dims < 1)
    throw Except("Grid not initialized", __FILE__, __LINE__);

  if (data_center_pos.size() != dims)
    throw Except("Improper data center position vector", __FILE__, __LINE__);

  if (data_spread.size() != dims)
    throw Except("Improper data spread vector", __FILE__, __LINE__);

  // CHECK SPREAD VECTOR ELEMENTS???

  // CHECK WEIGHT???

  // Get the total volume of the data region

  double data_vol = 1.0;

  for (unsigned int i = 0; i < dims; i++)
    data_vol *= data_spread[i];

  if (!(data_vol > 0.0))
    add_data(data_center_pos, data_value, weight);
  else
    {
    const double oneOver_data_vol = 1.0 / data_vol;

    // Determine the range of bins

    std::vector<double> min_data_pos = data_center_pos;
    std::vector<double> max_data_pos = data_center_pos;

    for (unsigned int i = 0; i < dims; i++)
      {
      const double axis_half_width = 0.5 * data_spread[i];

      min_data_pos[i] -= axis_half_width;
      max_data_pos[i] += axis_half_width;
      }

    std::vector<unsigned int> min_bin_index;
    std::vector<unsigned int> max_bin_index;

    grid_info.get_nearest_cell(min_data_pos, min_bin_index);
    grid_info.get_nearest_cell(max_data_pos, max_bin_index);

    std::vector<unsigned int> bin_index(dims);

    // Set up the counter

    for (unsigned int i = 0; i < dims; i++)
      bin_index[i] = (max_bin_index[i] - min_bin_index[i]) + 1;

    VectorCounter bin_index_counter(bin_index);

    // Run through the overlapping bins

    std::vector<double> bin_min_coords;
    std::vector<double> bin_max_coords;

    while (bin_index_counter.in_bounds())
      {
      for (unsigned int i = 0; i < dims; i++)
	bin_index[i] = min_bin_index[i] + bin_index_counter[i];

      grid_info.get_cell_bounds(bin_index, bin_min_coords, bin_max_coords);

      const double overlap_volume 
	= Block::Utils::get_volume_overlap(bin_min_coords, bin_max_coords,
					   min_data_pos, max_data_pos);

      const double bin_weight = (overlap_volume * oneOver_data_vol) * weight;

      grid_data[bin_index].add_entry(data_value, bin_weight);

      bin_index_counter.increment();
      }
    }
  }


bool BinGrid::in_bounds(const std::vector<double> & coordinate,
			const std::vector<double> & spread) const
  {
  bool ret_val = false;

  // ERROR CHECKING???

  if (dims > 0)
    {
    ret_val = true;

    const std::vector<double> & grid_min_bounds = get_min_boundary();
    const std::vector<double> & grid_max_bounds = get_max_boundary();

    for (unsigned int i = 0; i < dims; i++)
      {
      const double axis_coord = coordinate[i];

      const double axis_half_spread = 0.5 * spread[i];

      const double min_coord = axis_coord - axis_half_spread;
      const double max_coord = axis_coord + axis_half_spread;

      const double axis_grid_min_coord = grid_min_bounds[i];
      const double axis_grid_max_coord = grid_max_bounds[i];

      if ((!(max_coord > axis_grid_min_coord)) 
	  || (!(min_coord < axis_grid_max_coord)))
	ret_val = false;
      }
    }

  return ret_val;
  }


}  // End namespace Stats
}  // End namespace QuickFlash
