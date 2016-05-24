// C++ header file quickflash_stats_bins.hpp

/*
  By Nathan C. Hearn
     February 18, 2008

  Statistics routines for binned data.
*/


#ifndef QUICKFLASH_STATS_BINS_HPP
#define QUICKFLASH_STATS_BINS_HPP


#include <vector>
#include <cmath>
#include "quickflash_block_blockinfo.hpp"
#include "quickflash_block_blockdata.hpp"


namespace QuickFlash
{
namespace Stats
{

// Class Bin

class Bin
  {
  public :

    Bin() :
      num_entries(0), sum_weight(0.0), sum_weight_values(0.0),
      sum_weight_sqValues(0.0), min_value(0.0), max_value(0.0)
      { }

    Bin(const Bin & source) :
      num_entries(0), sum_weight(0.0), sum_weight_values(0.0),
      sum_weight_sqValues(0.0), min_value(0.0), max_value(0.0)
      { reset(source); }

    ~Bin() { }

    Bin & operator=(const Bin & source)
      {
      reset(source);
      return *this;
      }

    Bin & operator+=(const Bin & source)
      {
      add_bin(source);
      return *this;
      }

    void reset();

    void reset(const Bin & source);

    void add_bin(const Bin & source);

    void add_entry(const double value, const double weight=1.0);      

    unsigned int get_num_entries() const { return num_entries; }

    bool entry_present() const
      { 
      bool present = false;
      
      if (num_entries > 0)
	present = true;

      return present;
      }

    double get_weight_sum() const { return sum_weight; }

    double get_mean() const 
      { 
      double mean_val = 0.0;

      if (sum_weight > 0.0)
	mean_val = sum_weight_values / sum_weight;

      return mean_val;
      }

    double get_weight_value_sum() const { return sum_weight_values; }
    double get_weight_sqValue_sum() const { return sum_weight_sqValues; }

    double get_pop_variance() const;
    double get_pop_stddev() const { return std::sqrt(get_pop_variance()); }

    double get_pop_mean_square_difference(const double reference_value=0.0) 
      const;
    double get_pop_root_mean_square_difference(const double 
					         reference_value=0.0) const
      { return std::sqrt(get_pop_mean_square_difference(reference_value)); }

    double get_min_value() const { return min_value; }
    double get_max_value() const { return max_value; }

  private :

    unsigned int num_entries;

    double sum_weight;

    double sum_weight_values;
    double sum_weight_sqValues;

    double min_value;
    double max_value;
  };


// Class BinArray

class BinArray
  {
  public :

    BinArray() :
      bin_count(0), min_array_bounds(0.0), max_array_bounds(0.0), 
      bin_min_bounds(), bin_center_coords(), bin_widths(), bins()
      { }

    BinArray(const double min_boundary, const unsigned int num_bins,
	     const double bin_width) :
      bin_count(0), min_array_bounds(0.0), max_array_bounds(0.0), 
      bin_min_bounds(), bin_center_coords(), bin_widths(), bins()
      { reset(min_boundary, num_bins, bin_width); }

    BinArray(const double min_boundary, const double max_boundary, 
	     const unsigned int num_bins, const bool log_scale=false) :
      bin_count(0), min_array_bounds(0.0), max_array_bounds(0.0), 
      bin_min_bounds(), bin_center_coords(), bin_widths(), bins()
      { reset(min_boundary, max_boundary, num_bins, log_scale); }

    BinArray(const std::vector<double> bin_min_boundary_list, 
	     const double max_boundary) :
      bin_count(0), min_array_bounds(0.0), max_array_bounds(0.0), 
      bin_min_bounds(), bin_center_coords(), bin_widths(), bins()
      { reset(bin_min_boundary_list, max_boundary); }

    BinArray(const BinArray & source) :
      bin_count(0), min_array_bounds(0.0), max_array_bounds(0.0), 
      bin_min_bounds(), bin_center_coords(), bin_widths(), bins()
      { reset(source); }

    ~BinArray() { }

    BinArray & operator=(const BinArray & source)
      {
      reset(source);
      return *this;
      }

    const Bin & operator[](const unsigned int bin_index) const
      { return bins[bin_index]; }

    void reset();

    void reset(const double min_boundary, const unsigned int num_bins,
	       const double bin_width);

    void reset(const double min_boundary, const double max_boundary, 
	       const unsigned int num_bins, const bool log_scale=false);

    void reset(const std::vector<double> bin_min_boundary_list, 
	       const double max_boundary);

    void reset(const BinArray & source);

    void reset_bins();

    void add_data(const double coordinate, const double value,
		  const double weight=1.0);

    void add_data_spread(const double coordinate, const double value,
			 const double spread, const double weight=1.0);

    unsigned int get_num_bins() const { return bin_count; }

    bool in_bounds(const double coordinate, const double spread=0.0) const;

    unsigned int get_bin_index(const double coordinate) const;

    unsigned int get_nearest_bin_index(const double coordinate) const;

    bool in_bin(const unsigned int bin_index, const double coordinate, 
		const double spread=0.0) const;

    double get_overlap_area(const unsigned int bin_index, 
			    const double coordinate, const double spread) 
      const;

    const Bin & get_bin(const unsigned int bin_index) const 
      { return bins[bin_index]; }

    double get_min_boundary() const { return min_array_bounds; }
    double get_max_boundary() const { return max_array_bounds; }

    const std::vector<double> & get_bin_centers() const
      { return bin_center_coords; }

    double get_bin_center(const unsigned int bin_index) const
      { return bin_center_coords[bin_index]; }

    const std::vector<double> & get_bin_widths() const
      { return bin_widths; }

    double get_bin_width(const unsigned int bin_index) const
      { return bin_widths[bin_index]; }

    const std::vector<double> & get_bin_min_bounds() const
      { return bin_min_bounds; }

    double get_bin_min_bounds(const unsigned int bin_index) const
      { return bin_min_bounds[bin_index]; }

    double get_bin_max_bounds(const unsigned int bin_index) const
      { 
      double ret_val = max_array_bounds;

      if (bin_index < (bin_count - 1))
	ret_val = bin_min_bounds[bin_index + 1];

      return ret_val;
      }

    double get_peak_mean(unsigned int & bin_index) const;

    double get_peak_mean() const
      {
      unsigned int bin_index = 0;

      return get_peak_mean(bin_index);
      }

    double get_transition_pos(const double threshold_value,
			      const bool start_lower_bound=true) const;
				
    double get_first_pos_below(const double threshold_value,
			       const bool start_lower_bound=true) const;
				
    double get_first_pos_above(const double threshold_value,
			       const bool start_lower_bound=true) const;

  private :

    static double find_threshold_pos(const double lo_pos, const double hi_pos,
				     const double lo_val, const double hi_val,
				     const double threshold_value);

  private :

    unsigned int bin_count;

    double min_array_bounds;
    double max_array_bounds;

    std::vector<double> bin_min_bounds;

    std::vector<double> bin_center_coords;
    std::vector<double> bin_widths;

    std::vector<Bin> bins;
  };


// Class BinGrid

class BinGrid
  {
  public :

    BinGrid();

    BinGrid(const std::vector<double> & grid_min_coords, 
	    const std::vector<double> & bin_grid_spacing, 
	    const std::vector<unsigned int> & bin_grid_dims);

    BinGrid(const BinGrid & source);

    ~BinGrid() { }

    BinGrid & operator=(const BinGrid & source)
      {
      reset(source);
      return *this;
      }

    const Bin & operator[](const unsigned int bin_index) const
      { return get_bin(bin_index); }

    const Bin & operator[](const std::vector<unsigned int> & bin_index) const
      { return get_bin(bin_index); }

    const Bin & get_bin(const unsigned int bin_index) const
      { return grid_data[bin_index]; }

    const Bin & get_bin(const std::vector<unsigned int> & bin_index) const
      { return grid_data[bin_index]; }

    void reset();

    void reset(const std::vector<double> & grid_min_coords, 
	       const std::vector<double> & bin_grid_spacing, 
	       const std::vector<unsigned int> & bin_grid_dims);

    void reset(const BinGrid & source);

    void add_data(const std::vector<double> & data_pos, const double value,
		  const double weight=1.0)
      {
      if (in_bounds(data_pos))
	{
	const unsigned int bin_index = get_nearest_bin_index(data_pos);

	grid_data[bin_index].add_entry(value, weight);
	}
      }

    void add_data_spread(const std::vector<double> & data_center_pos,
			 const double data_value,
			 const std::vector<double> & data_spread,
			 const double weight=1.0);

    unsigned int get_num_bins() const { return grid_info.get_num_cells(); }

    const std::vector<unsigned int> & get_grid_dims() const
      { return grid_info.get_block_dims(); }

    bool in_bounds(const std::vector<double> & coordinate) const
      { return grid_info.in_block(coordinate); }

    bool in_bounds(const std::vector<double> & coordinate,
		   const std::vector<double> & spread) const;

    unsigned int get_nearest_bin_index(const std::vector<double> & coordinate)
      const
      { return grid_info.get_nearest_cell(coordinate); }

    bool in_bin(const unsigned int bin_index, 
		const std::vector<double> & coordinate) const
      { return grid_info.in_cell(bin_index, coordinate); }

    const std::vector<double> & get_min_boundary() const
      { return grid_info.get_min_coords(); }

    const std::vector<double> & get_max_boundary() const
      { return grid_info.get_max_coords(); }

  private :

    unsigned int dims;

    Block::BlockInfo grid_info;
    Block::BlockData<Bin> grid_data;
  };


}  // End namespace Stats
}  // End namespace QuickFlash


#endif  // QUICKFLASH_STATS_BINS_HPP
