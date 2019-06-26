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
// A simple program to demonstrate some basic shading techniques
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
#include <sstream>
#include <chrono>

#include "geometry.h"
#include "objects.h"
#include "lights.h"
#include "loader.h"
#include "utils.h"

struct Options
{
	uint32_t width = 1920;
	uint32_t height = 1080;
	float fov = 90.f;
	Vec3f backgroundColor = kDefaultBackgroundColor;
	Matrix44f cameraToWorld;
	float bias = 0.0001f;
	uint32_t maxDepth = 5;
};

struct IntersectInfo
{
	const Object *hitObject = nullptr;
	float tNear = kInfinity;
	Vec2f uv;
	uint32_t index = 0;
};

bool trace(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	IntersectInfo &isect,
	RayType rayType = kPrimaryRay)
{
	isect.hitObject = nullptr;
	for (uint32_t k = 0; k < objects.size(); ++k) {
		float tNear = kInfinity;
		uint32_t index = 0;
		Vec2f uv;
		if (objects[k]->intersect(orig, dir, tNear, index, uv) && tNear < isect.tNear) {
			isect.hitObject = objects[k].get();
			isect.tNear = tNear;
			isect.index = index;
			isect.uv = uv;
		}
	}

	return (isect.hitObject != nullptr);
}

Vec3f castRay(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	const std::vector<std::unique_ptr<Light>> &lights,
	const Options &options,
	const uint32_t & depth = 0)
{
	if (depth > options.maxDepth) return options.backgroundColor;
	Vec3f hitColor = 0;
	IntersectInfo isect;
	if (trace(orig, dir, objects, isect)) {
		// [comment]
		// Evaluate surface properties (P, N, texture coordinates, etc.)
		// [/comment]
		Vec3f hitPoint = orig + dir * isect.tNear;
		Vec3f hitNormal;
		Vec2f hitTexCoordinates;
		isect.hitObject->getSurfaceProperties(hitPoint, dir, isect.index, isect.uv, hitNormal, hitTexCoordinates);
		switch (isect.hitObject->materialType) {
		// [comment]
		// Simulate diffuse object
		// [/comment]
		case kPhong:
		{
			// [comment]
			// Light loop (loop over all lights in the scene and accumulate their contribution)
			// [/comment]
			Vec3f diffuse = 0, specular = 0;
			for (uint32_t i = 0; i < lights.size(); ++i) {
				Vec3f lightDir, lightIntensity;
				IntersectInfo isectShad;
				lights[i]->illuminate(hitPoint, lightDir, lightIntensity, isectShad.tNear);

				bool vis = !trace(hitPoint + hitNormal * options.bias, -lightDir, objects, isectShad, kShadowRay);

				// compute the diffuse component
				diffuse += vis * isect.hitObject->diffuseColor * lightIntensity * std::max(0.f, hitNormal.dotProduct(-lightDir));

				// compute the specular component
				// what would be the ideal reflection direction for this light ray
				Vec3f R = reflect(lightDir, hitNormal);
				specular += vis * lightIntensity * std::pow(std::max(0.f, R.dotProduct(-dir)), isect.hitObject->specularExponent);
			}
			hitColor = diffuse * isect.hitObject->Kd + specular * isect.hitObject->Ks;
			//std::cerr << hitColor << std::endl;
			break;
		}
		default:
			break;
		}
	}
	else {
		hitColor = options.backgroundColor;
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
	const std::vector<std::unique_ptr<Object>> &objects,
	const std::vector<std::unique_ptr<Light>> &lights)
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
			*(pix++) = castRay(orig, dir, objects, lights, options);
		}
		fprintf(stderr, "\r%3d%c", uint32_t(j / (float)options.height * 100), '%');
	}
	auto timeEnd = std::chrono::high_resolution_clock::now();
	auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
	fprintf(stderr, "\rDone: %.2f (sec)\n", passedTime / 1000);

	// save framebuffer to file
	std::ofstream ofs;
	ofs.open("ray_phong_shading.ppm");
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.height * options.width; ++i) {
		char r = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].x));
		char g = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].y));
		char b = (char)(255 * scratch::utils::clamp(0, 1, framebuffer[i].z));
		ofs << r << g << b;
	}
	ofs.close();
}

// [comment]
// In the main function of the program, we create the scene (create objects and lights)
// as well as set the options for the render (image widht and height, maximum recursion
// depth, field-of-view, etc.). We then call the render function().
// [/comment]
int main(int argc, char **argv)
{
	// loading gemetry
	std::vector<std::unique_ptr<Object>> objects;
	// lights
	std::vector<std::unique_ptr<Light>> lights;
	Options options;

	// aliasing example
	options.fov = 36.87f;
	options.width = 1920;
	options.height = 1080;
	options.cameraToWorld[3][2] = 12;
	options.cameraToWorld[3][1] = 1;

	Matrix44f xform;
	xform[0][0] = 1.f;
	xform[1][1] = 1.f;
	xform[2][2] = 1.f;
	TriangleMesh *mesh = scratch::loader::loadPolyMeshFromFile("data/plane.geo", &xform);
	if (mesh != nullptr) {
		mesh->smoothShading = false;
		objects.push_back(std::unique_ptr<Object>(mesh));
	}

	float w[5] = {0.04f, 0.08f, 0.1f, 0.15f, 0.2f};
	for (int i = -4, n = 2, k = 0; i <= 4; i+= 2, n *= 5, k++) {
		Matrix44f xformSphere;
		xformSphere[3][0] = i;
		xformSphere[3][1] = 1.f;
		Sphere *sph = new Sphere(xformSphere, 0.9);
		sph->diffuseColor = Vec3f(0.18f, 0.18f, 0.18f);
		sph->materialType = kPhong;
		sph->specularExponent = n;
		sph->Ks = w[k];
		objects.push_back(std::unique_ptr<Object>(sph));
	}

	Matrix44f l2w(11.146836f, -5.781569f, -0.0605886f, 0.f, -1.902827f, -3.543982f, -11.895445f, 0.f, 5.459804f, 10.568624f, -4.02205f, 0.f, 0.f, 0.f, 0.f, 1.f);
	lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1.f, 5.f)));

	// finally, render
	render(options, objects, lights);

	return 0;
}