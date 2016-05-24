// C++ header file quickflash_surface_isosurface.hpp

/*
  By Nathan C. Hearn
     July 16, 2008

  Class for encapsulating isosurface construction functions.
*/


#ifndef QUICKFLASH_SURFACE_ISOSURFACE_HPP
#define QUICKFLASH_SURFACE_ISOSURFACE_HPP


#include <list>
#include "quickflash_file_meshinfo.hpp"
#include "quickflash_file_dataset.hpp"
#include "quickflash_surface_triangle.hpp"

namespace QuickFlash
{
namespace Surface
{

// Class IsosurfaceGen

class IsosurfaceGen
  {
  public :

    IsosurfaceGen() : 
      meshinfo(0), dset(0), limit_refine_level(false), max_refine_level(0),
      limit_vol_min(false), vol_min_bounds(), limit_vol_max(false),
      vol_max_bounds()
      { }

    IsosurfaceGen(const File::MeshInfo & mesh_info,
		  const File::Dataset & isosurface_dataset,
		  const bool limit_refinement_level=false,
		  const unsigned int max_refinement_level=1) :
      meshinfo(0), dset(0), limit_refine_level(false), max_refine_level(0),
      limit_vol_min(false), vol_min_bounds(), limit_vol_max(false),
      vol_max_bounds()
      { 
      reset(mesh_info, isosurface_dataset, limit_refinement_level,
	    max_refinement_level);
      }

    IsosurfaceGen(const IsosurfaceGen & source) :
      meshinfo(0), dset(0), limit_refine_level(false), max_refine_level(0),
      limit_vol_min(false), vol_min_bounds(), limit_vol_max(false),
      vol_max_bounds()
      { reset(source); }

    ~IsosurfaceGen() { }

    IsosurfaceGen & operator=(const IsosurfaceGen & source)
      {
      reset(source);
      return *this;
      }

    void reset();

    void reset(const File::MeshInfo & mesh_info,
	       const File::Dataset & isosurface_dataset,
	       const bool limit_refinement_level=false,
	       const unsigned int max_refinement_level=1);

    void reset(const IsosurfaceGen & source);

    void set_refinement_limit(const unsigned int max_refinement_level=1,
			      const bool limit_refinement_level=true);

    void set_volume_limit(const bool limit_volume_min_bounds=false,
			  const std::vector<double> 
			    & volume_min_bounds=std::vector<double>(),
			  const bool limit_volume_max_bounds=false,
			  const std::vector<double> 
			  & volume_max_bounds=std::vector<double>());

    void add_block_surface(const unsigned int block_index,
			   const double isolevel_value,
			   std::list<Triangle> & triangle_list);

  private :

    bool block_in_bounds(const std::vector<double> & min_coords,
			 const std::vector<double> & max_coords) const;

    bool triangle_in_bounds(const Triangle & triangle) const;

  private :

    static const unsigned int dims;

    const File::MeshInfo * meshinfo;  // Do not delete
    const File::Dataset * dset;     // Do not delete

    bool limit_refine_level;
    unsigned int max_refine_level;

    bool limit_vol_min;
    std::vector<double> vol_min_bounds;

    bool limit_vol_max;
    std::vector<double> vol_max_bounds;
  };


}  // End namespace Surface
}  // End namespace QuickFlash


#endif  // QUICKFLASH_SURFACE_ISOSURFACE_HPP
