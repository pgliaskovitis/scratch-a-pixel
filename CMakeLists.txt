cmake_minimum_required (VERSION 3.1)
project (scratch_a_pixel)

# can change with cmake -DCMAKE_BUILD_TYPE=Debug ..
set (CMAKE_CXX_STANDARD 14)
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()

if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
	set (CMAKE_CXX_FLAGS "-Wall -Wextra -ffast-math")
	set (CMAKE_CXX_FLAGS_DEBUG "-g")
	set (CMAKE_CXX_FLAGS_RELEASE "-O3 -g")
elseif (${CMAKE_C_COMPILER_ID} STREQUAL "CLANG")
	set (CMAKE_CXX_FLAGS "-Weverything")
	set (CMAKE_CXX_FLAGS_DEBUG "-g")
	set (CMAKE_CXX_FLAGS_RELEASE "-O3")
elseif (${CMAKE_C_COMPILER_ID} STREQUAL "MSVC")
	set (CMAKE_CXX_FLAGS_RELEASE "/O2")
endif ()

include_directories (${CMAKE_SOURCE_DIR}/src/include)
if (WIN32)
	include_directories (${CMAKE_SOURCE_DIR})
endif ()

# basic color spaces/imaging
add_executable (img_mcbeth src/img/mcbeth.cpp)
add_executable (img_read_write src/img/read_write.cpp)
add_executable (img_bokeh src/img/bokeh.cpp)

# basic Monte-Carlo
add_executable (mc_samplingdist src/mc/sampling_dist.cpp)
if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
	add_executable (mc_inverse_transform src/mc/inverse_transform.cpp)
endif ()
add_executable (mc_sim src/mc/sim.cpp)
add_executable (mc_integration src/mc/integration.cpp)
add_executable (mc_random src/mc/random.cpp)
add_executable (mc_importance src/mc/importance.cpp)
add_executable (mc_quasi src/mc/quasi.cpp)

# basic rasterization
add_executable (geometry src/raster/geometry.cpp)
add_executable (persp_proj src/raster/persp_proj.cpp)
add_executable (pinhole src/raster/pinhole.cpp)
add_executable (gl_ortho_proj_matrix src/raster/gl_ortho_proj_matrix.cpp)
add_executable (gl_proj_matrix src/raster/gl_proj_matrix.cpp)
add_executable (proj_matrix src/raster/proj_matrix.cpp)
add_executable (raster_2d_triangle src/raster/2d_triangle.cpp)
add_executable (raster_3d_triangle src/raster/3d_triangle.cpp)
add_executable (raster_3d_mesh src/raster/3d_mesh.cpp)
add_executable (raster_3d_mesh_loader src/raster/3d_mesh_loader.cpp)

# basic raytracing
add_executable (ray_basic_tracer src/ray/basic_tracer.cpp)
add_executable (ray_whitted src/ray/whitted.cpp)
add_executable (ray_camera_rays src/ray/camera_rays.cpp)
add_executable (ray_simple_shapes src/ray/simple_shapes.cpp)
add_executable (ray_box src/ray/box.cpp)
add_executable (ray_triangle src/ray/triangle.cpp)
add_executable (ray_3d_mesh src/ray/3d_mesh.cpp)
add_executable (ray_3d_transform src/ray/3d_transform.cpp)
add_executable (ray_simple_shading src/ray/simple_shading.cpp)
add_executable (ray_phong_shading src/ray/phong_shading.cpp)
add_executable (ray_indirect_diffuse_shading src/ray/indirect_diffuse_shading.cpp)

#advanced raytracing
add_executable (ray_bezier_teapot src/ray/bezier_teapot.cpp)
add_executable (ray_acceleration src/ray/acceleration.cpp)
add_executable (ray_implicit_2D src/ray/2d_implicit.cpp)
add_executable (ray_distance_fields src/ray/distance_fields.cpp)

#basic volume rendering
add_executable (ray_skycolor src/ray/skycolor.cpp)

# basic noise
add_executable (noise_value src/noise/value.cpp)
add_executable (noise_perlin src/noise/perlin.cpp)
