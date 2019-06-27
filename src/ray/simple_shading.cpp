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
	float fov = 90;
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
			if (rayType == kShadowRay && objects[k]->materialType == kReflectionAndRefraction) continue;
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
			case kDiffuse:
			{
			// [comment]
			// Light loop (loop over all lights in the scene and accumulate their contribution)
			// [/comment]
			for (uint32_t i = 0; i < lights.size(); ++i) {
				Vec3f lightDir, lightIntensity;
				IntersectInfo isectShad;
				lights[i]->illuminate(hitPoint, lightDir, lightIntensity, isectShad.tNear);
				bool vis = !trace(hitPoint + hitNormal * options.bias, -lightDir, objects, isectShad, kShadowRay);
				// compute the pattern
				float angle = scratch::utils::deg2rad(45);
				float s = hitTexCoordinates.x * cos(angle) - hitTexCoordinates.y * sin(angle);
				float t = hitTexCoordinates.y * cos(angle) + hitTexCoordinates.x * sin(angle);
				float scaleS = 20, scaleT = 20;
				// float pattern = (cos(hitTexCoordinates.y * 2 * M_PI * scaleT) * sin(hitTexCoordinates.x * 2 * M_PI * scaleS) + 1) * 0.5; // isect.hitObject->albedo
				float pattern = (scratch::utils::modulo(s * scaleS) < 0.5) ^ (scratch::utils::modulo(t * scaleT) < 0.5);
				// float pattern = (scratch::utils::modulo(s * scaleS) < 0.5);
				hitColor += vis * pattern * lightIntensity * std::max(0.f, hitNormal.dotProduct(-lightDir));
			}
			break;
		}
		// [comment]
		// Simulate reflection only
		// [/comment]
		case kReflection:
		{
			Vec3f R = reflect(dir, hitNormal);
			R.normalize();
			break;
		}
		// [comment]
		// Simulate transparent object (reflection/transmission/fresnel)
		// [/comment]
		case kReflectionAndRefraction:
		{
			Vec3f refractionColor = 0, reflectionColor = 0;
			// compute fresnel
			float kr;
			fresnel(dir, hitNormal, isect.hitObject->ior, kr);
			bool outside = dir.dotProduct(hitNormal) < 0;
			Vec3f bias = options.bias * hitNormal;
			// compute refraction if it is not a case of total internal reflection
			if (kr < 1) {
				Vec3f refractionDirection = refract(dir, hitNormal, isect.hitObject->ior).normalize();
				Vec3f refractionRayOrig = outside ? hitPoint - bias : hitPoint + bias;
				refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, depth + 1);
			}

			Vec3f reflectionDirection = reflect(dir, hitNormal).normalize();
			Vec3f reflectionRayOrig = outside ? hitPoint + bias : hitPoint - bias;
			reflectionColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, depth + 1);

			// mix the two
			hitColor += reflectionColor * kr + refractionColor * (1 - kr);
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
	std::ofstream ofs;
	ofs.open("ray_simple_shading.ppm", std::ios::out | std::ios::binary);
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

	auto diffuseReflectionRefractionFunc = [&](void) {
		// glass and pen example
		// setting up options
		options.fov = 36.87f;
		options.maxDepth = 10;
		options.bias = 0.001f;
		options.width = 1920;
		options.height = 1080;
		options.cameraToWorld = Matrix44f(-0.972776f, 0.f, -0.231748f, 0.f, -0.114956f, 0.8683f, 0.482536f, 0.f, 0.201227f, 0.49604f, -0.844661f, 0.f, 6.696465f, 22.721296f, -30.097976f, 1.f);

		TriangleMesh *mesh1 = scratch::loader::loadPolyMeshFromFile("data/backdrop.geo", &Matrix44f::kIdentity);
		if (mesh1 != nullptr) {
			mesh1->materialType = kDiffuse;
			objects.push_back(std::unique_ptr<Object>(mesh1));
		}

		TriangleMesh *mesh3 = scratch::loader::loadPolyMeshFromFile("data/cylinder.geo", &Matrix44f::kIdentity);
		if (mesh3 != nullptr) {
			mesh3->materialType = kReflectionAndRefraction;
			mesh3->ior = 1.5f;
			objects.push_back(std::unique_ptr<Object>(mesh3));
		}

		TriangleMesh *mesh4 = scratch::loader::loadPolyMeshFromFile("data/pen.geo", &Matrix44f::kIdentity);
		if (mesh4 != nullptr) {
			mesh4->materialType = kDiffuse;
			mesh4->diffuseColor = Vec3f(0.18f, 0.18f, 0.18f);
			mesh4->smoothShading = false;
			objects.push_back(std::unique_ptr<Object>(mesh4));
		}

		Matrix44f xform1;
		xform1[3][0] = -1.2f;
		xform1[3][1] = 6.f;
		xform1[3][2] = -3.f;
		Sphere *sph1 = new Sphere(xform1, 5.f);
		sph1->materialType = kReflectionAndRefraction;

		Matrix44f l2w(11.146836f, -5.781569f, -0.0605886f, 0.f, -1.902827f, -3.543982f, -11.895445f, 0.f, 5.459804f, 10.568624f, -4.02205f, 0.f, 0.f, 0.f, 0.f, 1.f);
		lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 1)));

		// finally, render
		render(options, objects, lights);
	};

	auto patternsFunc = [&](void) {
		// simple plane example (patterns)
		options.fov = 36.87f;
		options.width = 1920;
		options.height = 1080;
		options.cameraToWorld = Matrix44f(0.707107f, 0.f, -0.707107f, 0.f, -0.331295f, 0.883452f, -0.331295f, 0.f, 0.624695f, 0.468521f, 0.624695f, 0.f, 28.f, 21.f, 28.f, 1.f);

		TriangleMesh *mesh = scratch::loader::loadPolyMeshFromFile("data/plane.geo", &Matrix44f::kIdentity);
		if (mesh != nullptr) {
			mesh->materialType = kDiffuse;
			mesh->diffuseColor = Vec3f(0.18f, 0.18f, 0.18f);
			mesh->smoothShading = false;
			objects.push_back(std::unique_ptr<Object>(mesh));
		}

		Matrix44f l2w(11.146836f, -5.781569f, -0.0605886f, 0.f, -1.902827f, -3.543982f, -11.895445f, 0.f, 5.459804f, 10.568624f, -4.02205f, 0.f, 0.f, 0.f, 0.f, 1.f);
		lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 1)));

		// finally, render
		render(options, objects, lights);
	};

	auto multipleReflectionRefractionFunc = [&](void) {
		// multiple glasses example
		options.fov = 36.87f;
		options.width = 1920;
		options.height = 1080;
		options.cameraToWorld = Matrix44f(0.999945f, 0.f, 0.0104718f, 0.f, 0.00104703f, 0.994989f, -0.0999803f, 0.f, -0.0104193f, 0.0999858f, 0.994934f, 0.f, -0.978596f, 17.911879f, 75.483369f, 1.f);

		TriangleMesh *mesh = scratch::loader::loadPolyMeshFromFile("data/glasses.geo", &Matrix44f::kIdentity);
		if (mesh != nullptr) {
			mesh->materialType = kReflectionAndRefraction;
			mesh->ior = 1.3f;
			objects.push_back(std::unique_ptr<Object>(mesh));
		}

		TriangleMesh *mesh1 = scratch::loader::loadPolyMeshFromFile("data/backdrop1.geo", &Matrix44f::kIdentity);
		if (mesh1 != nullptr) {
			mesh1->materialType = kDiffuse;
			mesh1->diffuseColor = Vec3f(0.18f, 0.18f, 0.18f);
			objects.push_back(std::unique_ptr<Object>(mesh1));
		}

		Matrix44f l2w(0.95292f, 0.289503f, 0.0901785f, 0.f, -0.0960954f, 0.5704f, -0.815727f, 0.f, -0.287593f, 0.768656f, 0.571365f, 0.f, 0.f, 0.f, 0.f, 1.f);
		lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 1)));

		// finally, render
		render(options, objects, lights);
	};

	auto aliasingFunc = [&](void) {
		// aliasing example
		options.fov = 36.87f;
		options.width = 1920;
		options.height = 1080;
		options.cameraToWorld = Matrix44f(0.999945f, 0.f, 0.0104718f, 0.f, 0.00104703f, 0.994989f, -0.0999803f, 0.f, -0.0104193f, 0.0999858f, 0.994934f, 0.f, -0.978596f, 17.911879f, 75.483369f, 1.f);

		Matrix44f xform;
		xform[0][0] = 10;
		xform[1][1] = 10;
		xform[2][2] = 10;
		xform[3][2] = -40;
		TriangleMesh *mesh = scratch::loader::loadPolyMeshFromFile("data/plane.geo", &xform);
		if (mesh != nullptr) {
			mesh->materialType = kDiffuse;
			mesh->diffuseColor = Vec3f(0.18f, 0.18f, 0.18f);
			mesh->smoothShading = false;
			objects.push_back(std::unique_ptr<Object>(mesh));
		}

		Matrix44f l2w(11.146836f, -5.781569f, -0.0605886f, 0.f, -1.902827f, -3.543982f, -11.895445f, 0.f, 5.459804f, 10.568624f, -4.02205f, 0.f, 0.f, 0.f, 0.f, 1.f);
		lights.push_back(std::unique_ptr<Light>(new DistantLight(l2w, 1, 1)));

		// finally, render
		render(options, objects, lights);
	};

#if 0
	diffuseReflectionRefractionFunc();
#elif 0
	patternsFunc();
#elif 1
	multipleReflectionRefractionFunc();
#else
	aliasingFunc();
#endif

	return 0;
}
