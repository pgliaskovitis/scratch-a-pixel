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

# basic Monte-Carlo
add_executable (sampling_dist sampling_dist.cpp)
if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
	add_executable (inversetransform inversetransform.cpp)
endif ()
add_executable (mcsim mcsim.cpp)
add_executable (mcintegration mcintegration.cpp)
add_executable (mcrandom mcrandom.cpp)
add_executable (mcimportance mcimportance.cpp)
add_executable (mcquasi mcquasi.cpp)

# basic raytracing
add_executable (raytracer raytracer.cpp)

# basic rasterization
add_executable (perspproj perspproj.cpp)
add_executable (pinhole pinhole.cpp)
add_executable (raster2d_triangle raster2d_triangle.cpp)
add_executable (raster3d_triangle raster3d_triangle.cpp)
add_executable (raster3d raster3d.cpp)
add_executable (glorthoprojmatrix glorthoprojmatrix.cpp)
add_executable (glprojmatrix glprojmatrix.cpp)
add_executable (projmatrix projmatrix.cpp)

# basic color spaces/imaging
add_executable (mcbeth mcbeth.cpp)
add_executable (readwrite_img readwrite_img.cpp)
add_executable (bokeh bokeh.cpp)
