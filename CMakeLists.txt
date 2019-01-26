cmake_minimum_required (VERSION 3.1)
project (scratch_a_pixel)

set (CMAKE_CXX_STANDARD 11)
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Debug)
endif ()

if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
	set (CMAKE_CXX_FLAGS "-Wall -Wextra")
	set (CMAKE_CXX_FLAGS_DEBUG "-g")
	set (CMAKE_CXX_FLAGS_RELEASE "-O3")
elseif (${CMAKE_C_COMPILER_ID} STREQUAL "MSVC")
	set (CMAKE_CXX_FLAGS_RELEASE "/Ox")
endif ()

include_directories (${CMAKE_SOURCE_DIR}/src/include)
if (WIN32)
	include_directories (${CMAKE_SOURCE_DIR})
endif ()

# basic color spaces/imaging
add_executable (img_mcbeth src/img/mcbeth.cpp)
add_executable (img_readwrite src/img/readwrite.cpp)
add_executable (img_bokeh src/img/bokeh.cpp)

# basic Monte-Carlo
add_executable (mc_samplingdist src/mc/samplingdist.cpp)
if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
	add_executable (mc_inversetransform src/mc/inversetransform.cpp)
endif ()
add_executable (mc_sim src/mc/sim.cpp)
add_executable (mc_integration src/mc/integration.cpp)
add_executable (mc_random src/mc/random.cpp)
add_executable (mc_importance src/mc/importance.cpp)
add_executable (mc_quasi src/mc/quasi.cpp)

# basic rasterization
add_executable (geometry src/raster/geometry.cpp)
add_executable (perspproj src/raster/perspproj.cpp)
add_executable (pinhole src/raster/pinhole.cpp)
add_executable (glorthoprojmatrix src/raster/glorthoprojmatrix.cpp)
add_executable (glprojmatrix src/raster/glprojmatrix.cpp)
add_executable (projmatrix src/raster/projmatrix.cpp)
add_executable (raster2d_triangle src/raster/raster2d_triangle.cpp)
add_executable (raster3d_triangle src/raster/raster3d_triangle.cpp)
add_executable (raster3d src/raster/raster3d.cpp)

# basic raytracing
add_executable (ray_basictracer src/ray/basictracer.cpp)
add_executable (ray_whitted src/ray/whitted.cpp)
add_executable (ray_camerarays src/ray/camerarays.cpp)
add_executable (ray_simpleshapes src/ray/simpleshapes.cpp)
add_executable (ray_box src/ray/box.cpp)
add_executable (ray_triangle src/ray/triangle.cpp)
