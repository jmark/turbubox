include ../settings.make


OBJ = quickflash_file_flashdefs.o quickflash_file_dataset.o \
	quickflash_file_datafile.o quickflash_hdf5.o \
	quickflash_block_blockinfo.o \
	quickflash_file_siminfo.o quickflash_file_meshinfo.o \
	quickflash_file_treeinfo.o \
	quickflash_block_blocktype.o quickflash_tree_builder.o \
	quickflash_block_utils.o quickflash_tree_block.o \
	quickflash_counters.o \
	quickflash_geometry.o quickflash_block_cellvolume.o \
	quickflash_particles_data.o quickflash_mesh_boundarydefs.o \
	quickflash_tree_nodedata.o quickflash_tree_nodeid.o \
	quickflash_slice_slicegrid.o quickflash_slice_slicer.o \
	quickflash_slice_output.o \
	quickflash_color_color.o quickflash_color_colormap.o \
	quickflash_color_utils.o \
	quickflash_image_imagegrid.o quickflash_image_imagemagick.o \
	quickflash_utils_text.o quickflash_utils_commandline.o \
	quickflash_utils_file.o quickflash_stats_bins.o \
        quickflash_surface_triangle.o \
	quickflash_surface_marchingtetrahedrons.o \
	quickflash_surface_isosurface.o


LIB_FILES = libquickflash.a libquickflash.$(SO_EXT)


CCFLAGS += $(CCFLAGS_SHARED)  # Does this interfere with libquickflash.a???


all : $(LIB_FILES)

clean : clean-obj
	rm -f $(LIB_FILES)

clean-obj :
	rm -f *.o

libquickflash.a : $(OBJ)
	$(AR) $(ARFLAGS) libquickflash.a $(OBJ)
	$(RANLIB) libquickflash.a

libquickflash.$(SO_EXT) : $(OBJ)
	$(LD) -o libquickflash.$(SO_EXT) $(OBJ) $(LDFLAGS) $(LDFLAGS_SHARED) \
	    $(LDFLAGS_HDF5) $(LDFLAGS_MAGICK)

quickflash_file_flashdefs.o : quickflash_file_flashdefs.hpp \
	    quickflash_file_flashdefs.cpp quickflash_hdf5.hpp \
	    quickflash_except.hpp 
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_file_flashdefs.o \
	    quickflash_file_flashdefs.cpp

quickflash_file_dataset.o : quickflash_file_dataset.hpp \
	    quickflash_file_dataset.cpp \
	    quickflash_block_blockdata.hpp quickflash_block_blockdata.tcc \
	    quickflash_except.hpp quickflash_counters.hpp \
	    quickflash_cache_indexcache.hpp quickflash_cache_indexcache.tcc \
	    quickflash_cache_item.hpp quickflash_hdf5.hpp \
	    quickflash_file_flashdefs.hpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_file_dataset.o \
	    quickflash_file_dataset.cpp

quickflash_file_datafile.o : quickflash_file_datafile.hpp \
	    quickflash_file_datafile.cpp quickflash_file_siminfo.hpp \
	    quickflash_file_flashdefs.hpp quickflash_geometry.hpp \
	    quickflash_file_meshinfo.hpp quickflash_block_utils.hpp \
	    quickflash_block_blocktype.hpp quickflash_tree_block.hpp \
	    quickflash_tree_node.hpp quickflash_tree_node.tcc \
	    quickflash_tree_nodeid.hpp quickflash_tree_nodedata.hpp \
	    quickflash_block_blockinfo.hpp quickflash_block_cellvolume.hpp \
	    quickflash_block_domain.hpp quickflash_except.hpp \
	    quickflash_mesh_boundarydefs.hpp quickflash_file_dataset.hpp \
	    quickflash_block_blockdata.hpp quickflash_block_blockdata.tcc \
	    quickflash_counters.hpp quickflash_cache_indexcache.hpp \
	    quickflash_cache_indexcache.tcc quickflash_cache_item.hpp \
	    quickflash_registry.hpp quickflash_particles_data.hpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_file_datafile.o \
	    quickflash_file_datafile.cpp


quickflash_hdf5.o : quickflash_hdf5.hpp quickflash_hdf5.cpp 
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_hdf5.o \
	    quickflash_hdf5.cpp

quickflash_counters.o : quickflash_counters.hpp quickflash_counters.cpp
	$(CXX) -c $(CCFLAGS) -o quickflash_counters.o quickflash_counters.cpp

quickflash_block_blockinfo.o : quickflash_block_blockinfo.hpp \
	    quickflash_block_blockinfo.cpp \
	    quickflash_except.hpp quickflash_block_utils.hpp \
	    quickflash_geometry.hpp quickflash_block_cellvolume.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_block_blockinfo.o \
	    quickflash_block_blockinfo.cpp

quickflash_block_blocktype.o : quickflash_block_blocktype.hpp \
	    quickflash_block_blocktype.cpp
	$(CXX) -c $(CCFLAGS) -o quickflash_block_blocktype.o \
	    quickflash_block_blocktype.cpp

quickflash_block_utils.o : quickflash_block_utils.hpp \
	    quickflash_block_utils.cpp quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_block_utils.o -c quickflash_block_utils.cpp

quickflash_tree_block.o : quickflash_tree_block.hpp quickflash_tree_block.cpp \
	    quickflash_tree_node.hpp quickflash_tree_node.tcc \
	    quickflash_tree_nodeid.hpp \
	    quickflash_tree_nodedata.hpp quickflash_block_blockinfo.hpp \
	    quickflash_block_utils.hpp quickflash_block_domain.hpp \
	    quickflash_geometry.hpp quickflash_block_cellvolume.hpp \
	    quickflash_mesh_boundarydefs.hpp quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_tree_block.o \
	    quickflash_tree_block.cpp

quickflash_tree_builder.o : quickflash_tree_builder.hpp \
	    quickflash_tree_builder.cpp \
	    quickflash_tree_block.hpp quickflash_tree_node.hpp \
	    quickflash_tree_node.tcc quickflash_tree_nodeid.hpp \
	    quickflash_tree_nodedata.hpp \
	    quickflash_block_blockinfo.hpp quickflash_block_utils.hpp \
	    quickflash_block_domain.hpp quickflash_except.hpp \
	    quickflash_geometry.hpp quickflash_block_cellvolume.hpp \
	    quickflash_mesh_boundarydefs.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_tree_builder.o \
	    quickflash_tree_builder.cpp

quickflash_tree_nodeid.o : quickflash_tree_nodeid.hpp \
	    quickflash_tree_nodeid.cpp
	$(CXX) -c $(CCFLAGS) -o quickflash_tree_nodeid.o \
	    quickflash_tree_nodeid.cpp

quickflash_file_treeinfo.o : quickflash_file_treeinfo.hpp \
	    quickflash_file_treeinfo.cpp quickflash_except.hpp \
	    quickflash_hdf5.hpp \
	    quickflash_mesh_boundarydefs.hpp quickflash_file_flashdefs.hpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_file_treeinfo.o \
	    quickflash_file_treeinfo.cpp

quickflash_file_meshinfo.o : quickflash_file_meshinfo.hpp \
	    quickflash_file_meshinfo.cpp quickflash_file_siminfo.hpp \
	    quickflash_file_flashdefs.hpp quickflash_geometry.hpp \
	    quickflash_block_utils.hpp quickflash_block_blocktype.hpp \
	    quickflash_tree_block.hpp quickflash_tree_node.hpp \
	    quickflash_tree_node.tcc quickflash_tree_nodeid.hpp \
	    quickflash_tree_nodedata.hpp quickflash_block_blockinfo.hpp \
	    quickflash_block_cellvolume.hpp quickflash_block_domain.hpp \
	    quickflash_except.hpp quickflash_mesh_boundarydefs.hpp \
	    quickflash_file_treeinfo.hpp quickflash_tree_builder.hpp 
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_file_meshinfo.o \
	    quickflash_file_meshinfo.cpp

quickflash_file_siminfo.o : quickflash_file_siminfo.hpp \
	    quickflash_file_siminfo.cpp \
	    quickflash_file_flashdefs.hpp quickflash_geometry.hpp \
	    quickflash_except.hpp quickflash_hdf5.hpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_file_siminfo.o \
	    quickflash_file_siminfo.cpp

quickflash_geometry.o : quickflash_geometry.hpp quickflash_geometry.cpp
	$(CXX) -c $(CCFLAGS) -o quickflash_geometry.o quickflash_geometry.cpp

quickflash_block_cellvolume.o : quickflash_block_cellvolume.hpp \
	    quickflash_block_cellvolume.cpp quickflash_geometry.hpp \
	    quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_block_cellvolume.o \
	    quickflash_block_cellvolume.cpp

quickflash_particles_data.o : quickflash_particles_data.hpp \
	    quickflash_particles_data.cpp quickflash_hdf5.hpp \
	    quickflash_file_flashdefs.hpp quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_particles_data.o \
	    quickflash_particles_data.cpp

quickflash_mesh_boundarydefs.o : quickflash_mesh_boundarydefs.hpp \
	    quickflash_mesh_boundarydefs.cpp
	$(CXX) -c $(CCFLAGS) -o quickflash_mesh_boundarydefs.o \
	    quickflash_mesh_boundarydefs.cpp

quickflash_tree_nodedata.o : quickflash_tree_nodedata.hpp \
	    quickflash_tree_nodedata.cpp quickflash_block_blockinfo.hpp \
	    quickflash_block_utils.hpp quickflash_geometry.hpp \
	    quickflash_block_cellvolume.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_tree_nodedata.o \
	    quickflash_tree_nodedata.cpp

quickflash_slice_slicegrid.o : quickflash_slice_slicegrid.hpp \
	    quickflash_slice_slicegrid.cpp quickflash_block_blockinfo.hpp \
	    quickflash_block_blockdata.hpp quickflash_block_blockdata.tcc \
	    quickflash_block_utils.hpp \
	    quickflash_geometry.hpp quickflash_block_cellvolume.hpp \
	    quickflash_except.hpp quickflash_counters.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_slice_slicegrid.o \
	    quickflash_slice_slicegrid.cpp

quickflash_slice_slicer.o : quickflash_slice_slicer.hpp \
	    quickflash_slice_slicer.cpp \
	    quickflash_file_meshinfo.hpp quickflash_file_siminfo.hpp \
	    quickflash_file_flashdefs.hpp quickflash_geometry.hpp \
	    quickflash_block_utils.hpp quickflash_block_blocktype.hpp \
	    quickflash_tree_block.hpp quickflash_tree_node.hpp \
	    quickflash_tree_node.tcc quickflash_tree_nodeid.hpp \
	    quickflash_tree_nodedata.hpp quickflash_block_blockinfo.hpp \
	    quickflash_block_cellvolume.hpp quickflash_block_domain.hpp \
	    quickflash_except.hpp quickflash_geometry.hpp \
	    quickflash_mesh_boundarydefs.hpp quickflash_file_dataset.hpp \
	    quickflash_block_blockdata.hpp quickflash_block_blockdata.tcc \
	    quickflash_counters.hpp quickflash_cache_indexcache.hpp \
	    quickflash_cache_indexcache.tcc quickflash_cache_item.hpp \
	    quickflash_slice_slicegrid.hpp 
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) -o quickflash_slice_slicer.o \
	    quickflash_slice_slicer.cpp

quickflash_slice_output.o : quickflash_slice_output.hpp \
	    quickflash_slice_output.cpp \
	    quickflash_slice_slicegrid.hpp quickflash_block_blockinfo.hpp \
	    quickflash_block_blockdata.hpp quickflash_block_blockdata.tcc \
	    quickflash_block_utils.hpp quickflash_geometry.hpp \
	    quickflash_block_cellvolume.hpp quickflash_except.hpp \
	    quickflash_counters.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_slice_output.o \
	    quickflash_slice_output.cpp

quickflash_color_color.o : quickflash_color_color.hpp quickflash_color_color.cpp \
	    quickflash_color_utils.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_color_color.o \
	    quickflash_color_color.cpp

quickflash_color_colormap.o : quickflash_color_colormap.hpp \
	    quickflash_color_colormap.cpp quickflash_color_color.hpp \
	    quickflash_except.hpp quickflash_utils_text.hpp \
	    quickflash_utils_text.tcc
	$(CXX) -c $(CCFLAGS) -o quickflash_color_colormap.o \
	    quickflash_color_colormap.cpp

quickflash_color_utils.o : quickflash_color_utils.hpp \
	    quickflash_color_utils.cpp
	$(CXX) -c $(CCFLAGS) -o quickflash_color_utils.o \
	    quickflash_color_utils.cpp

quickflash_image_imagegrid.o : quickflash_image_imagegrid.hpp \
	    quickflash_image_imagegrid.cpp quickflash_block_blockdata.hpp \
	    quickflash_block_blockdata.tcc quickflash_except.hpp \
	    quickflash_counters.hpp quickflash_color_color.hpp \
	    quickflash_color_colormap.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_image_imagegrid.o \
	    quickflash_image_imagegrid.cpp

quickflash_image_imagemagick.o : quickflash_image_imagemagick.hpp \
	    quickflash_image_imagemagick.cpp \
	    quickflash_image_imagegrid.hpp quickflash_except.hpp \
	    quickflash_color_color.hpp quickflash_block_blockdata.hpp \
	    quickflash_block_blockdata.tcc quickflash_color_colormap.hpp \
	  quickflash_counters.hpp
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_MAGICK) \
	    -o quickflash_image_imagemagick.o quickflash_image_imagemagick.cpp

quickflash_utils_text.o : quickflash_utils_text.hpp quickflash_utils_text.tcc \
	    quickflash_utils_text.cpp quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_utils_text.o \
	    quickflash_utils_text.cpp

quickflash_utils_commandline.o : quickflash_utils_commandline.hpp \
	    quickflash_utils_commandline.cpp quickflash_utils_option.hpp \
	    quickflash_utils_text.hpp quickflash_utils_text.tcc \
	    quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_utils_commandline.o \
	    quickflash_utils_commandline.cpp

quickflash_utils_file.o : quickflash_utils_file.hpp quickflash_utils_file.cpp
	$(CXX) -c $(CCFLAGS) -o quickflash_utils_file.o \
	    quickflash_utils_file.cpp

quickflash_stats_bins.o : quickflash_stats_bins.hpp quickflash_stats_bins.cpp \
	    quickflash_block_blockinfo.hpp quickflash_block_utils.hpp \
	    quickflash_geometry.hpp quickflash_block_cellvolume.hpp \
	    quickflash_block_blockdata.hpp quickflash_block_blockdata.tcc \
	    quickflash_counters.hpp quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_stats_bins.o \
	    quickflash_stats_bins.cpp

quickflash_surface_triangle.o : quickflash_surface_triangle.hpp \
	    quickflash_surface_triangle.cpp quickflash_geometry_vector.hpp \
	    quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_surface_triangle.o \
	    quickflash_surface_triangle.cpp

quickflash_surface_marchingtetrahedrons.o : \
	    quickflash_surface_marchingtetrahedrons.hpp \
	    quickflash_surface_marchingtetrahedrons.cpp \
	    quickflash_surface_triangle.hpp quickflash_geometry_vector.hpp \
	    quickflash_surface_cube.hpp \
	    quickflash_surface_tetrahedron.hpp quickflash_except.hpp
	$(CXX) -c $(CCFLAGS) -o quickflash_surface_marchingtetrahedrons.o \
	    quickflash_surface_marchingtetrahedrons.cpp

quickflash_surface_isosurface.o : quickflash_surface_isosurface.hpp \
	    quickflash_surface_isosurface.cpp quickflash_file_meshinfo.hpp \
	    quickflash_file_siminfo.hpp quickflash_file_flashdefs.hpp \
	    quickflash_geometry.hpp quickflash_block_utils.hpp \
	    quickflash_block_blocktype.hpp quickflash_tree_block.hpp \
	    quickflash_tree_node.hpp quickflash_tree_node.tcc \
	    quickflash_tree_nodeid.hpp quickflash_tree_nodedata.hpp \
	    quickflash_block_blockinfo.hpp quickflash_block_cellvolume.hpp \
	    quickflash_block_domain.hpp quickflash_except.hpp \
	    quickflash_mesh_boundarydefs.hpp quickflash_file_dataset.hpp \
	    quickflash_block_blockdata.hpp quickflash_block_blockdata.tcc \
	    quickflash_counters.hpp quickflash_cache_indexcache.hpp \
	    quickflash_cache_indexcache.tcc quickflash_cache_item.hpp \
	    quickflash_surface_triangle.hpp quickflash_geometry_vector.hpp \
	    quickflash_surface_cube.hpp \
	    quickflash_surface_marchingtetrahedrons.hpp \
	    quickflash_surface_tetrahedron.hpp 
	$(CXX) -c $(CCFLAGS) $(CCFLAGS_HDF5) \
	    -o quickflash_surface_isosurface.o \
	    quickflash_surface_isosurface.cpp
