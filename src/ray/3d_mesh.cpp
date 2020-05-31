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
// A simple program to demonstrate how to ray-trace a polygon mesh
//[/header]

#include <fstream>
#include <chrono>

#include "objects.h"
#include "generator.h"
#include "loader.h"
#include "utils.h"

struct Options
{
	uint32_t width = 1920;
	uint32_t height = 1080;
	float fov = 90.f;
	Vec3f backgroundColor = kDefaultBackgroundColor;
	Matrix44f cameraToWorld;
};

bool trace(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	float &tNear, uint32_t &index, Vec2f &uv, Object **hitObject)
{
	*hitObject = nullptr;
	for (uint32_t k = 0; k < objects.size(); ++k) {
		float tNearTriangle = kInfinity;
		uint32_t indexTriangle;
		Vec2f uvTriangle;
		if (objects[k]->intersect(orig, dir, tNearTriangle, indexTriangle, uvTriangle) && tNearTriangle < tNear) {
			*hitObject = objects[k].get();
			tNear = tNearTriangle;
			index = indexTriangle;
			uv = uvTriangle;
		}
	}

	return (*hitObject != nullptr);
}

Vec3f castRay(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	const Options &options)
{
	Vec3f hitColor = options.backgroundColor;
	float tnear = kInfinity;
	Vec2f uv;
	uint32_t index = 0;
	Object *hitObject = nullptr;
	if (trace(orig, dir, objects, tnear, index, uv, &hitObject)) {
		Vec3f hitPoint = orig + dir * tnear;
		Vec3f hitNormal;
		Vec2f hitTexCoordinates;
		hitObject->getSurfaceProperties(hitPoint, dir, index, uv, hitNormal, hitTexCoordinates);
		float NdotView = std::max(0.f, hitNormal.dotProduct(-dir));
		const int M = 10;
		float checker = (fmod(hitTexCoordinates.x * M, 1.0f) > 0.5f) ^ (fmod(hitTexCoordinates.y * M, 1.0) < 0.5);
		float c = 0.3f * (1.f - checker) + 0.7f * checker;

		hitColor = c * NdotView; //Vec3f(uv.x, uv.y, 0);
	}

	return hitColor;
}

//[comment]
// The main render function. This where we iterate over all pixels in the image, generate primary rays and cast these rays
// into the scene. The content of the framebuffer is saved to a file.
//[/comment]
void render(
	const Options &options,
	const std::vector<std::unique_ptr<Object>> &objects,
	const uint32_t &frame)
{
	std::unique_ptr<Vec3f []> framebuffer(new Vec3f[options.width * options.height]);
	Vec3f *pix = framebuffer.get();
	float scale = tan(scratch::utils::deg2rad(options.fov * 0.5f));
	float imageAspectRatio = options.width / (float)options.height;
	Vec3f orig;
	options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
	auto timeStart = std::chrono::high_resolution_clock::now();
	for (uint32_t j = 0; j < options.height; ++j) {
		for (uint32_t i = 0; i < options.width; ++i) {
			// generate primary ray direction
			float x = (2 * (i + 0.5f) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5f) / (float)options.height) * scale;
			Vec3f dir;
			options.cameraToWorld.multDirMatrix(Vec3f(x, y, -1), dir);
			dir.normalize();
			*(pix++) = castRay(orig, dir, objects, options);
		}
		fprintf(stderr, "\r%3d%c", uint32_t(j / (float)options.height * 100), '%');
	}
	auto timeEnd = std::chrono::high_resolution_clock::now();
	auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
	std::cerr << std::endl << "Wall passed time:  " << passedTime << " ms" << std::endl;

	// save framebuffer to file
	char buff[256];
	sprintf(buff, "ray_3d_mesh.%04d.ppm", frame);
	std::ofstream ofs;
	ofs.open(buff, std::ios::out | std::ios::binary);
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.height * options.width; ++i) {
		char r = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].x));
		char g = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].y));
		char b = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].z));
		ofs << r << g << b;
	}
	ofs.close();
}

//[comment]
// In the main function of the program, we create the scene (create objects and lights) as well as set the options for the
// render (image widht and height, maximum recursion depth, field-of-view, etc.). We then call the render function().
//[/comment]
int main(int argc, char **argv)
{
	// setting up options
	Options options;
	//options.cameraToWorld[3][2] = 10;
	Matrix44f tmp = Matrix44f(0.707107f, -0.331295f, 0.624695f, 0.f, 0.f, 0.883452f, 0.468521f, 0.f, -0.707107f, -0.331295f, 0.624695f, 0.f, -1.63871f, -5.747777f, -40.400412f, 1.f);
	options.cameraToWorld = tmp.inverse();
	options.fov = 50.0393f;
#if 1
	std::vector<std::unique_ptr<Object>> objects;
	TriangleMesh *mesh = scratch::loader::loadPolyMeshFromFile("data/cow.geo", nullptr);
	if (mesh != nullptr) objects.push_back(std::unique_ptr<Object>(mesh));

	// finally, render
	render(options, objects, 0);
#else
	for (uint32_t i = 0; i < 10; ++i) {
		int divs = 5 + i;
		// creating the scene (adding objects and lights)
		std::vector<std::unique_ptr<Object>> objects;
		TriangleMesh *mesh = scratch::generator::generatePolySphere(11, divs);
		objects.push_back(std::unique_ptr<Object>(mesh));
		auto timeStart = std::chrono::high_resolution_clock::now();
		// finally, render
		render(options, objects, i);
		auto timeEnd = std::chrono::high_resolution_clock::now();
		auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
		std::cerr << mesh->numTris << " " << passedTime << std::endl;
	}
#endif

	return 0;
}