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
// This program simulates indirect diffuse using Monte-Carlo integration
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

#include <random>

#include "geometry.h"
#include "objects.h"
#include "lights.h"
#include "loader.h"
#include "utils.h"

static std::default_random_engine generator;
static std::uniform_real_distribution<float> distribution(0, 1);

struct Options
{
	uint32_t width = 1920;
	uint32_t height = 1080;
	float fov = 90.f;
	Vec3f backgroundColor = kDefaultBackgroundColor;
	Matrix44f cameraToWorld;
	float bias = 0.0001f;
	uint32_t maxDepth = 2;
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

void createCoordinateSystem(const Vec3f &N, Vec3f &Nt, Vec3f &Nb)
{
	if (std::fabs(N.x) > std::fabs(N.y)) {
		Nt = Vec3f(N.z, 0.f, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
	} else {
		Nt = Vec3f(0.f, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
	}
	Nb = N.crossProduct(Nt);
}

Vec3f uniformSampleHemisphere(const float &r1, const float &r2)
{
	// cos(theta) = u1 = y
	// cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
	float sinTheta = sqrtf(1.f - r1 * r1);
	float phi = 2 * M_PI * r2;
	float x = sinTheta * cosf(phi);
	float z = sinTheta * sinf(phi);
	return Vec3f(x, r1, z);
}

Vec3f castRay(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	const std::vector<std::unique_ptr<Light>> &lights,
	const Options &options,
	const uint32_t & depth = 0)
{
	if (depth > options.maxDepth) return 0;//options.backgroundColor;
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
		case kDiffuse:
		{
			// [comment]
			// Compute direct ligthing
			// [/comment]
			Vec3f directLighting = 0;
			for (uint32_t i = 0; i < lights.size(); ++i) {
				Vec3f lightDir, lightIntensity;
				IntersectInfo isectShad;
				lights[i]->illuminate(hitPoint, lightDir, lightIntensity, isectShad.tNear);
				bool vis = !trace(hitPoint + hitNormal * options.bias, -lightDir, objects, isectShad, kShadowRay);
				directLighting = vis * lightIntensity * std::max(0.f, hitNormal.dotProduct(-lightDir));
			}

			// [comment]
			// Compute indirect ligthing
			// [/comment]
			Vec3f indirectLigthing = 0;

			uint32_t N = 128;// / (depth + 1);
			Vec3f Nt, Nb;
			createCoordinateSystem(hitNormal, Nt, Nb);
			float pdf = 1 / (2 * M_PI);
			for (uint32_t n = 0; n < N; ++n) {
				float r1 = distribution(generator);
				float r2 = distribution(generator);
				Vec3f sample = uniformSampleHemisphere(r1, r2);
				Vec3f sampleWorld(
					sample.x * Nb.x + sample.y * hitNormal.x + sample.z * Nt.x,
					sample.x * Nb.y + sample.y * hitNormal.y + sample.z * Nt.y,
					sample.x * Nb.z + sample.y * hitNormal.z + sample.z * Nt.z);
				// don't forget to divide by PDF and multiply by cos(theta)
				indirectLigthing += r1 * castRay(hitPoint + sampleWorld * options.bias,
				sampleWorld, objects, lights, options, depth + 1) / pdf;
			}
			// divide by N
			indirectLigthing /= (float)N;

			hitColor = (directLighting / M_PI + 2 * indirectLigthing) * isect.hitObject->diffuseColor;
			break;
		}
		default:
			break;
		}
	} else {
		hitColor = 1;
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
	float scale = tan(scratch::utils::deg2rad(options.fov * 0.5));
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
	float gamma = 1;
	std::ofstream ofs;
	ofs.open("ray_indirect_diffuse_shading.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.height * options.width; ++i) {
		char r = (char)(255 * scratch::utils::clamp(0, 1, powf(framebuffer[i].x, 1/gamma)));
		char g = (char)(255 * scratch::utils::clamp(0, 1, powf(framebuffer[i].y, 1/gamma)));
		char b = (char)(255 * scratch::utils::clamp(0, 1, powf(framebuffer[i].z, 1/gamma)));
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
	options.fov = 39.89f;
	options.width = 1920;
	options.height = 1080;
	options.cameraToWorld = Matrix44f(0.965926f, 0.f, -0.258819f, 0.f, 0.0066019f, 0.999675f, 0.0246386f, 0.f, 0.258735f, -0.0255078f, 0.965612f, 0.f, 0.764985f, 0.791882f, 5.868275f, 1.f);

	TriangleMesh *plane = scratch::loader::loadPolyMeshFromFile("data/planegi.geo", &Matrix44f::kIdentity);
	if (plane != nullptr) {
		plane->diffuseColor = Vec3f(0.225f, 0.144f, 0.144f);
		plane->materialType = kDiffuse;
		plane->specularExponent = 10;
		plane->smoothShading = true;
		objects.push_back(std::unique_ptr<Object>(plane));
	}

	TriangleMesh *cube = scratch::loader::loadPolyMeshFromFile("data/cubegi.geo", &Matrix44f::kIdentity);
	if (cube != nullptr) {
		cube->diffuseColor = Vec3f(0.188559f, 0.287f, 0.200726f);
		cube->materialType = kDiffuse;
		cube->specularExponent = 10;
		cube->smoothShading = true;
		objects.push_back(std::unique_ptr<Object>(cube));
	}

	Matrix44f xformSphere;
	xformSphere[3][1] = 1;
	Sphere *sph = new Sphere(xformSphere, 1);
	if (sph != nullptr) {
		sph->diffuseColor = Vec3f(0.18f, 0.18f, 0.18f);
		sph->materialType = kDiffuse;
		sph->specularExponent = 10;
		objects.push_back(std::unique_ptr<Object>(sph));
	}

	Matrix44f l2w(0.916445f, -0.218118f, 0.335488f, 0.f, 0.204618f, -0.465058f, -0.861309f, 0.f, 0.343889f, 0.857989f, -0.381569f, 0.f, 0.f, 0.f, 0.f, 1.f);
	// lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 16)));

	// finally, render
	render(options, objects, lights);

	return 0;
}

