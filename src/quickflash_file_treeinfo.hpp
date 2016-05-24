// C++ header file quickflash_file_treeinfo.hpp

/*
  By Nathan C. Hearn
     October 4, 2008

  Reader for block centers, widths, and topology from Flash data files.
*/


#ifndef QUICKFLASH_FILE_TREEINFO_HPP
#define QUICKFLASH_FILE_TREEINFO_HPP


#include <hdf5.h>
#include <vector>
#include "quickflash_mesh_boundarydefs.hpp"
#include "quickflash_file_flashdefs.hpp"


namespace QuickFlash
{
namespace File
{

class TreeInfo
  {
  public :

    typedef std::vector<double> Vector_Double;
    typedef std::vector<Vector_Double> Vector_Vector_Double;

    typedef std::vector<unsigned int> Vector_UInt;
    typedef std::vector<Vector_UInt> Vector_Vector_UInt;

    typedef std::vector<Mesh::NeighborType> Vector_NeighborType;
    typedef std::vector<Vector_NeighborType> Vector_Vector_NeighborType;

    typedef std::vector<Mesh::BoundaryType> Vector_BoundaryType;
    typedef std::vector<Vector_BoundaryType> Vector_Vector_BoundaryType;

  public :

    TreeInfo();

    TreeInfo(const hid_t file_id, const FlashVersion flash_version, 
	     const unsigned int mesh_dims);

    TreeInfo(const TreeInfo & source);

    ~TreeInfo() { }

    TreeInfo & operator=(const TreeInfo & source)
      {
      reset(source);
      return *this;
      }

    void reset();
    void reset(const hid_t file_id, const FlashVersion flash_version, 
	       const unsigned int mesh_dims);
    void reset(const TreeInfo & source);

    unsigned int get_dims() const { return dims; }

    unsigned int get_num_blocks() const { return num_blocks; }

    const Vector_UInt & get_process_id() const { return block_process_id; }

    const Vector_Vector_Double & get_block_center() const
      { return block_center; }

    const Vector_Vector_Double & get_block_width() const
      { return block_width; }

    const Vector_UInt & get_block_parent_index() const
      { return block_parent_index; }

    const Vector_Vector_UInt & get_block_neighbor_indexes_lo() const
      { return block_neighbor_indexes_lo; }

    const Vector_Vector_UInt & get_block_neighbor_indexes_hi() const
      { return block_neighbor_indexes_hi; }

    const Vector_Vector_NeighborType & get_block_neighbor_types_lo() const
      { return block_neighbor_types_lo; }

    const Vector_Vector_NeighborType & get_block_neighbor_types_hi() const
      { return block_neighbor_types_hi; }

    const Vector_Vector_BoundaryType & get_block_boundary_types_lo() const
      { return block_boundary_types_lo; }

    const Vector_Vector_BoundaryType & get_block_boundary_types_hi() const
      { return block_boundary_types_hi; }

    const Vector_Vector_UInt & get_block_child_indexes() const
      { return block_child_indexes; }

    const Vector_UInt & get_block_type() const
      { return block_type; }

    const Vector_UInt & get_block_refine_level() const
      { return block_refine_level; }

    unsigned int get_min_refine_level() const { return min_refine_level; }

    unsigned int get_max_refine_level() const { return max_refine_level; }

  private :

    void read_block_info(const hid_t file_id, const unsigned int mesh_dims); 
    void read_neighbor_info(const FlashVersion flash_version,
			    const hid_t file_id);
    void read_refinement_levels(const hid_t file_id);
    void read_block_types(const hid_t file_id);
    void read_process_ids(const hid_t file_id);

  private :

    unsigned int dims;

    unsigned int num_blocks;

    std::vector< std::vector<double> > block_center;
    std::vector< std::vector<double> > block_width;

    std::vector<unsigned int> block_parent_index;

    std::vector< std::vector<unsigned int> > block_neighbor_indexes_lo;
    std::vector< std::vector<unsigned int> > block_neighbor_indexes_hi;

    std::vector< std::vector<Mesh::NeighborType> > block_neighbor_types_lo;
    std::vector< std::vector<Mesh::NeighborType> > block_neighbor_types_hi;

    std::vector< std::vector<Mesh::BoundaryType> > block_boundary_types_lo;
    std::vector< std::vector<Mesh::BoundaryType> > block_boundary_types_hi;

    std::vector< std::vector<unsigned int> > block_child_indexes;

    std::vector<unsigned int> block_type;

    std::vector<unsigned int> block_refine_level;

    unsigned int min_refine_level;
    unsigned int max_refine_level;

    std::vector<unsigned int> block_process_id;
  };


}  // End namespace File
}  // End namespace QuickFlash


#endif  // QUICKFLASH_FILE_TREEINFO_HPP
