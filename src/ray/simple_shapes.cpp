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
// A simple program that uses ray-tracing to render a scene made out of spheres
//[/header]

#include <vector>
#include <fstream>
#include <chrono>

#include "objects.h"
#include "utils.h"

#define MAYA_STYLE 0

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0, 1);

struct Options
{
	uint32_t width;
	uint32_t height;
	float fov;
	Matrix44f cameraToWorld;
};

// [comment]
// Returns true if the ray intersects an object. The variable tNear is set to the closest intersection distance and hitObject
// is a pointer to the intersected object. The variable tNear is set to infinity and hitObject is set null if no intersection
// was found.
// [/comment]
bool trace(const Vec3f &orig, const Vec3f &dir, const std::vector<std::unique_ptr<Object>> &objects, float &tNear, const Object *&hitObject)
{
	tNear = kInfinity;
	std::vector<std::unique_ptr<Object>>::const_iterator iter = objects.begin();
	for (; iter != objects.end(); ++iter) {
		float t = kInfinity;
		uint32_t indexK;
		Vec2f uvK;
		if ((*iter)->intersect(orig, dir, t, indexK, uvK) && t < tNear) {
			hitObject = iter->get();
			tNear = t;
		}
	}

	return (hitObject != nullptr);
}

// [comment]
// Compute the color at the intersection point if any (returns background color otherwise)
// [/comment]
Vec3f castRay(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects)
{
	Vec3f hitColor = 0;
	const Object *hitObject = nullptr; // this is a pointer to the hit object
	float t; // this is the intersection distance from the ray origin to the hit point
	if (trace(orig, dir, objects, t, hitObject)) {
		Vec3f Phit = orig + dir * t;
		Vec3f In;
		Vec3f Nhit;
		Vec2f uv;
		uint32_t index = 0;
		Vec2f tex;
		hitObject->getSurfaceProperties(Phit, In, index, uv, Nhit, tex);
		// Use the normal and texture coordinates to shade the hit point.
		// The normal is used to compute a simple facing ratio and the texture coordinate
		// to compute a basic checker board pattern
		float scale = 4;
		float pattern = (fmodf(tex.x * scale, 1.f) > 0.5f) ^ (fmodf(tex.y * scale, 1.f) > 0.5f);
		hitColor = std::max(0.f, Nhit.dotProduct(-dir)) * mix(hitObject->diffuseColor, hitObject->diffuseColor * 0.8f, pattern);
	}

	return hitColor;
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void render(
	const Options &options,
	const std::vector<std::unique_ptr<Object>> &objects)
{
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
	options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
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
			options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
			dir.normalize();
			*(pix++) = castRay(orig, dir, objects);
		}
	}

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream ofs("./ray_simpleshapes.ppm", std::ios::out | std::ios::binary);
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
// In the main function of the program, we create the scene (create objects)
// as well as set the options for the render (image widht and height etc.).
// We then call the render function().
// [/comment]
int main(int argc, char **argv)
{
	// creating the scene (adding objects and lights)
	std::vector<std::unique_ptr<Object>> objects;

	// generate a scene made of random spheres
	uint32_t numSpheres = 32;
	gen.seed(331);
	for (uint32_t i = 0; i < numSpheres; ++i) {
		Vec3f randPos((0.5f - dis(gen)) * 10.f, (0.5f - dis(gen)) * 10.f, (0.5f + dis(gen) * 10.f));
		float randRadius = (0.5f + dis(gen) * 0.5f);
		objects.push_back(std::unique_ptr<Object>(new Sphere(randPos, randRadius)));
	}

	// setting up options
	Options options;
	options.width = 1920;
	options.height = 1080;
	options.fov = 51.52f;
	options.cameraToWorld = Matrix44f(0.945519f, 0.f, -0.325569f, 0.f, -0.179534f, 0.834209f, -0.521403f, 0.f, 0.271593f, 0.551447f, 0.78876f, 0.f, 4.208271f, 8.374532f, 17.932925f, 1.f);

	auto t_start = std::chrono::high_resolution_clock::now();

	// finally, render
	render(options, objects);

	auto t_end = std::chrono::high_resolution_clock::now();
	auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	std::cerr << "Wall passed time:  " << passedTime << " ms" << std::endl;

	return 0;
}
