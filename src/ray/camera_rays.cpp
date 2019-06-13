/*
 * Copyright (C) 2012  www.scratchapixel.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//[header]
// A simple program to demonstrate how to implement Whitted-style ray-tracing
//[/header]

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include "geometry.h"
#include "utils.h"

#define MAYA_STYLE 0

struct Options
{
	uint32_t width;
	uint32_t height;
	float fov;
};

class Light
{
public:
	Light() {}
};

class Object
{
 public:
	Object() {}
	virtual ~Object() {}
};

// [comment]
// This function doesn't do much at the moment. It simply takes the ray direction
// and turn it into a color. Ray direction coordinates are un the range [-1,1].
// To normalized them, we just add 1 and divide the result by 2.
// [/comment]
Vec3f castRay(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	const std::vector<std::unique_ptr<Light>> &lights,
	const Options &options,
	uint32_t depth)
{
	Vec3f hitColor = (dir + Vec3f(1.f)) * 0.5f;
	return hitColor;
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void render(
	const Options &options,
	const std::vector<std::unique_ptr<Object>> &objects,
	const std::vector<std::unique_ptr<Light>> &lights)
{
	Matrix44f cameraToWorld;
	Vec3f *framebuffer = new Vec3f[options.width * options.height];
	Vec3f *pix = framebuffer;
	float scale = tan(scratch::utils::deg2rad(options.fov * 0.5f));
	float imageAspectRatio = options.width / (float)options.height;
	// [comment]
	// Don't forget to transform the ray origin (which is also the camera origin
	// by transforming the point with coordinates (0,0,0) to world-space using the
	// camera-to-world matrix.
	// [/comment]
	Vec3f orig;
	cameraToWorld.multVecMatrix(Vec3f(0), orig);
	for (uint32_t j = 0; j < options.height; ++j) {
		for (uint32_t i = 0; i < options.width; ++i) {
			// [comment]
			// Generate primary ray direction. Compute the x and y position
			// of the ray in screen space. This gives a point on the image plane
			// at z=1. From there, we simply compute the direction by normalized
			// the resulting vec3f variable. This is similar to taking the vector
			// between the point on the image plane and the camera origin, which
			// in camera space is (0,0,0):
			//
			// ray.dir = normalize(Vec3f(x,y,-1) - Vec3f(0));
			// [/comment]
#ifdef MAYA_STYLE
			float x = (2 * (i + 0.5f) / (float)options.width - 1) * scale;
			float y = (1 - 2 * (j + 0.5f) / (float)options.height) * scale * 1 / imageAspectRatio;
#elif
			float x = (2 * (i + 0.5f) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5f) / (float)options.height) * scale;
#endif
			// [comment]
			// Don't forget to transform the ray direction using the camera-to-world matrix.
			// [/comment]
			Vec3f dir;
			cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
			dir.normalize();
			*(pix++) = castRay(orig, dir, objects, lights, options, 0);
		}
	}

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream ofs("./ray_camerarays.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.height * options.width; ++i) {
		char r = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].x));
		char g = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].y));
		char b = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].z));
		ofs << r << g << b;
	}

	ofs.close();

	delete [] framebuffer;
}

// [comment]
// In the main function of the program, we create the scene (create objects and lights)
// as well as set the options for the render (image widht and height, maximum recursion
// depth, field-of-view, etc.). We then call the render function().
// [/comment]
int main(int argc, char **argv)
{
	// creating the scene (adding objects and lights)
	std::vector<std::unique_ptr<Object>> objects;
	std::vector<std::unique_ptr<Light>> lights;

	// setting up options
	Options options;
	options.width = 1920;
	options.height = 1080;
	options.fov = 90.f;

	// finally, render
	render(options, objects, lights);

	return 0;
}
