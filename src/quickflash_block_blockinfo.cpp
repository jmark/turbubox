// C++ program file quickflash_block_blockinfo.cpp

/*
  By Nathan C. Hearn
     September 2, 2006

  Meta-info for block.
*/


#include "quickflash_block_blockinfo.hpp"
#include <vector>
#include "quickflash_geometry.hpp"
#include "quickflash_except.hpp"
#include "quickflash_block_cellvolume.hpp"


namespace QuickFlash
{
namespace Block
{

// Class BlockInfo


void BlockInfo::reset()
  {
  dims = 0;

  numcells = 0;

  blockdims.clear();

  mincoords.clear();
  maxcoords.clear();

  cellwidth.clear();
  oneOver_cellwidth.clear();
  cellvolume = 0.0;

  geom_type = Geometry::Cartesian;

  if (phys_cell_vol != 0)
    {
    delete phys_cell_vol;
    phys_cell_vol = 0;
    }
  }

    
void BlockInfo::reset(const unsigned int space_dims, 
		      const double block_mincoords[],
		      const double block_maxcoords[], 
		      const unsigned int block_dims[],
		      const Geometry::GeometryType geometry_type)
  {
  reset();  // Deletes phys_cell_vol

  dims = space_dims;

  numcells = 1;  // For starters...

  blockdims.resize(dims);

  mincoords.resize(dims);
  maxcoords.resize(dims);

  cellwidth.resize(dims);
  oneOver_cellwidth.resize(dims);

  for (unsigned int i = 0; i < dims; i++)
    {
    const unsigned int axis_cells = block_dims[i];

    if (axis_cells < 1)
      throw Except("Incompatible axis cell count", __FILE__, __LINE__);

    numcells *= axis_cells;

    const double min_axis_coord = block_mincoords[i];
    const double max_axis_coord = block_maxcoords[i];

    const double axis_width = max_axis_coord - min_axis_coord;

    if (!(axis_width > 0.0))
      throw Except("Incompatible block width", __FILE__, __LINE__);
  
    const double axis_cell_width 
      = axis_width / static_cast<double>(axis_cells);

    blockdims[i] = axis_cells;

    mincoords[i] = min_axis_coord;
    maxcoords[i] = max_axis_coord;

    cellwidth[i] = axis_cell_width;
    oneOver_cellwidth[i] = 1.0 / axis_cell_width;
    }

  geom_type = geometry_type;

  compute_cell_volume();
  }


void BlockInfo::reset(const std::vector<double> & block_mincoords,
		      const std::vector<double> & block_maxcoords,
		      const std::vector<unsigned int> & block_dims,
		      const Geometry::GeometryType geometry_type)
  {
  reset();  // Deletes phys_cell_vol

  dims = block_dims.size();

  if (block_mincoords.size() != dims)
    throw Except("Incompatible block min coords vector", __FILE__, __LINE__);

  if (block_maxcoords.size() != dims)
    throw Except("Incompatible block max coords vector", __FILE__, __LINE__);

  numcells = 1;  // For starters...

  blockdims.resize(dims);

  mincoords.resize(dims);
  maxcoords.resize(dims);

  cellwidth.resize(dims);
  oneOver_cellwidth.resize(dims);

  for (unsigned int i = 0; i < dims; i++)
    {
    const unsigned int axis_cells = block_dims[i];

    if (axis_cells < 1)
      throw Except("Incompatible axis cell count", __FILE__, __LINE__);

    numcells *= axis_cells;

    const double min_axis_coord = block_mincoords[i];
    const double max_axis_coord = block_maxcoords[i];

    const double axis_width = max_axis_coord - min_axis_coord;

    if (!(axis_width > 0.0))
      throw Except("Incompatible block width", __FILE__, __LINE__);
  
    const double axis_cell_width 
      = axis_width / static_cast<double>(axis_cells);

    blockdims[i] = axis_cells;

    mincoords[i] = min_axis_coord;
    maxcoords[i] = max_axis_coord;

    cellwidth[i] = axis_cell_width;
    oneOver_cellwidth[i] = 1.0 / axis_cell_width;
    }

  geom_type = geometry_type;

  compute_cell_volume();
  }


void BlockInfo::reset(const BlockInfo & source)
  {
  if (&source != this)
    {
    reset();  // Deletes phys_cell_vol

    dims = source.dims;
    numcells = source.numcells;
    blockdims = source.blockdims;

    mincoords = source.mincoords;
    maxcoords = source.maxcoords;

    cellwidth = source.cellwidth;
    oneOver_cellwidth = source.oneOver_cellwidth;
    cellvolume = source.cellvolume;

    geom_type = source.geom_type;

    if (source.phys_cell_vol != 0)
      phys_cell_vol = (source.phys_cell_vol)->clone();
    }
  }


void BlockInfo::get_cell_bounds(const std::vector<unsigned int> & cell_index,
				std::vector<double> & min_coords,
				std::vector<double> & max_coords) const
  {
  min_coords = mincoords;
  max_coords = mincoords;

  for (unsigned int i = 0; i < dims; i++)
    {
    // Determine the minimum coordinates for this cell

    const unsigned int axis_index = cell_index[i];

    // Keep boundaries from cell to cell consistent

    double min_cell_coord = 0.0;
    double max_cell_coord = 0.0;

    get_axis_cell_bounds(i, axis_index, min_cell_coord, max_cell_coord);

    // Make the upper bound of the last cell same as the upper block bound

    const unsigned int axis_next_index = axis_index + 1;

    if (axis_next_index >= blockdims[i])
      max_cell_coord = maxcoords[i];

    min_coords[i] = min_cell_coord;
    max_coords[i] = max_cell_coord;
    }
  }


unsigned int BlockInfo::get_nearest_cell(const std::vector<double> & coords)
  const
  {
  unsigned int cell_index = 0;

  // Construct the index axis-by-axis -- row-major storage

  if (dims > 0)
    {
    cell_index = get_nearest_cell_axis_index(0, coords[0]);

    for (unsigned int i = 1; i < dims; i++)
      {
      cell_index *= blockdims[i];
      cell_index += get_nearest_cell_axis_index(i, coords[i]);
      }
    }

  return cell_index;
  }


void BlockInfo::get_nearest_cell(const std::vector<double> & coords,
				 std::vector<unsigned int> & cell_index) const
  {
  cell_index.resize(dims);

  for (unsigned int i = 0; i < dims; i++)
    cell_index[i] = get_nearest_cell_axis_index(i, coords[i]);
  }


void BlockInfo::get_vector_index(const unsigned int scalar_index,
				 std::vector<unsigned int> & vector_index) 
  const
  {
  // No error checking??

  vector_index.resize(dims);

  if (dims > 0)
    {
    unsigned int work_index = scalar_index;

    for (unsigned int i = (dims - 1); i > 0; i--)
      {
      const unsigned int axis_blockdims = blockdims[i];

      vector_index[i] = work_index % axis_blockdims;
      work_index /= axis_blockdims;
      }

    vector_index[0] = work_index;
    }
  }


unsigned int BlockInfo::get_scalar_index(const std::vector<unsigned int> 
					 & vector_index) const
  {
  // No error checking??

  unsigned int index = 0;

  // Construct the index using row-major ordering

  if (dims > 0)
    {
    index = vector_index[0];

    for (unsigned int i = 1; i < dims; i++)
      {
      index *= blockdims[i];
      index += vector_index[i];
      }
    }

  return index;
  }


unsigned int BlockInfo::get_nearest_cell_axis_index(const unsigned int axis, 
						    const double axis_pos) 
  const
  {
  /*
    Find the axis cell index of the nearest cell to axis_pos.

    Note: The cell boundaries are defined with get_axis_cell_bounds
    for consistency.
  */

  const unsigned int block_axis_size = blockdims[axis];

  const double block_min_bound = mincoords[axis];

  const double rel_axis_coord = axis_pos - block_min_bound;

  const double oneOver_width = oneOver_cellwidth[axis];

  unsigned int axis_index = 0;

  if (!(rel_axis_coord < 0.0))
    {
    axis_index = static_cast<unsigned int>(rel_axis_coord * oneOver_width);

    if (axis_index >= block_axis_size)
      axis_index = block_axis_size - 1;
    }

  // Make sure this index is compatible with the cell boundaries

  double cell_min_bound = 0.0;
  double cell_max_bound = 0.0;

  get_axis_cell_bounds(axis, axis_index, cell_min_bound, cell_max_bound);

  if (axis_index < (block_axis_size - 1))
    {
    // Check upper-bound overlap (consistent with get_axis_cell_bounds)

    if (!(cell_max_bound > axis_pos))
      axis_index++;
    }

  if (axis_index > 0)
    {
    // Check lower-bound overlap (consistent with get_axis_cell_bounds)

    if (cell_min_bound > axis_pos)
      axis_index--;
    }

  return axis_index;
  }


void BlockInfo::compute_cell_volume()
  {
  // WARNING: All other object data members must be set up before calling!!!

  // Compute the coordinate volume

  double vol = 1.0;

  for (unsigned int i = 0; i < dims; i++)
    vol *= cellwidth[i];

  cellvolume = vol;

  // Build the physical volume object

  std::vector<double> min_cell_center(dims);

  for (unsigned int i = 0; i < dims; i++)
    min_cell_center[i] = mincoords[i] + (0.5 * cellwidth[i]);

  // Just in case one exists, delete the CellVolume object

  if (phys_cell_vol != 0)
    {
    delete phys_cell_vol;
    phys_cell_vol = 0;
    }

  phys_cell_vol = CellVolume::new_volume_obj(geom_type, blockdims,
					     cellwidth, min_cell_center);

  if (phys_cell_vol == 0)
    throw Except("Unable to build CellVolume object", __FILE__, __LINE__);
  }


void BlockInfo::get_axis_cell_bounds(const unsigned int axis,
				     const unsigned int axis_cell_index,
				     double & min_bound, double & max_bound) 
  const
  {
  /*
    Find cell boundaries along a given axis.

    Note: These boundaries are meant to be consistent from cell to cell
  */

  const double block_min_bound = mincoords[axis];
  const double width = cellwidth[axis];

  min_bound = block_min_bound 
    + (width * static_cast<double>(axis_cell_index));

  if (axis_cell_index < (blockdims[axis] - 1))
    max_bound = block_min_bound 
      + (width * static_cast<double>(axis_cell_index + 1));
  else
    max_bound = maxcoords[axis];
  }


}  // End namespace Block
}  // End namespace QuickFlash
