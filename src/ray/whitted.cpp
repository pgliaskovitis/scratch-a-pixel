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

#include <vector>
#include <fstream>
#include <chrono>

#include "objects.h"
#include "lights.h"

struct Options
{
	uint32_t width;
	uint32_t height;
	float fov;
	float imageAspectRatio;
	uint8_t maxDepth;
	Vec3f backgroundColor;
	float bias;
};

struct State
{
	uint32_t numPrimaryRays;
	uint32_t numReflectionRays;
	uint32_t numRefractionRays;
	uint32_t numShadowRays;
};

// [comment]
// Returns true if the ray intersects an object, false otherwise.
//
// \param orig is the ray origin
//
// \param dir is the ray direction
//
// \param objects is the list of objects the scene contains
//
// \param[out] tNear contains the distance to the cloesest intersected object.
//
// \param[out] index stores the index of the intersect triangle if the interesected object is a mesh.
//
// \param[out] uv stores the u and v barycentric coordinates of the intersected point
//
// \param[out] *hitObject stores the pointer to the intersected object (used to retrieve material information, etc.)
//
// \param isShadowRay is it a shadow ray. We can return from the function sooner as soon as we have found a hit.
// [/comment]
bool trace(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	float &tNear, uint32_t &index, Vec2f &uv, Object **hitObject)
{
	*hitObject = nullptr;
	for (uint32_t k = 0; k < objects.size(); ++k) {
		float tNearK = kInfinity;
		uint32_t indexK;
		Vec2f uvK;
		if (objects[k]->intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear) {
			*hitObject = objects[k].get();
			tNear = tNearK;
			index = indexK;
			uv = uvK;
		}
	}

	return (*hitObject != nullptr);
}

// [comment]
// Implementation of the Whitted-syle light transport algorithm (E [S*] (D|G) L)
//
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
//
// If the material of the intersected object is either reflective or reflective and refractive,
// then we compute the reflection/refracton direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refraction depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is diffuse/glossy we use the Phong illumation model to compute the color
// at the intersection point.
// [/comment]
Vec3f castRay(
	const Vec3f &orig, const Vec3f &dir,
	const std::vector<std::unique_ptr<Object>> &objects,
	const std::vector<std::unique_ptr<LightLite>> &lights,
	const Options &options,
	State &state,
	uint32_t depth,
	bool test = false)
{
	if (depth > options.maxDepth) {
		return options.backgroundColor;
	}

	Vec3f hitColor = options.backgroundColor;
	float tnear = kInfinity;
	Vec2f uv;
	uint32_t index = 0;
	Object *hitObject = nullptr;
	if (trace(orig, dir, objects, tnear, index, uv, &hitObject)) {
		Vec3f hitPoint = orig + dir * tnear;
		Vec3f N; // normal
		Vec2f st; // st coordinates
		hitObject->getSurfaceProperties(hitPoint, dir, index, uv, N, st);
		Vec3f tmp = hitPoint;
		switch (hitObject->materialType) {
			case kReflectionAndRefraction:
			{
				Vec3f reflectionDirection = reflect(dir, N).normalize();
				Vec3f refractionDirection = refract(dir, N, hitObject->ior).normalize();
				Vec3f reflectionRayOrig = (reflectionDirection.dotProduct(N) < 0) ?
					hitPoint - N * options.bias :
					hitPoint + N * options.bias;
				Vec3f refractionRayOrig = (refractionDirection.dotProduct(N) < 0) ?
					hitPoint - N * options.bias :
					hitPoint + N * options.bias;
				state.numReflectionRays++;
				Vec3f reflectionColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, state, depth + 1, 1);
				state.numRefractionRays++;
				Vec3f refractionColor = castRay(refractionRayOrig, refractionDirection, objects, lights, options, state, depth + 1, 1);
				float kr;
				fresnel(dir, N, hitObject->ior, kr);
				hitColor = reflectionColor * kr + refractionColor * (1 - kr);
				break;
			}
			case kReflection:
			{
				float kr;
				fresnel(dir, N, hitObject->ior, kr);
				Vec3f reflectionDirection = reflect(dir, N);
				Vec3f reflectionRayOrig = (reflectionDirection.dotProduct(N) < 0) ?
					hitPoint + N * options.bias :
					hitPoint - N * options.bias;
				state.numReflectionRays++;
				hitColor = castRay(reflectionRayOrig, reflectionDirection, objects, lights, options, state, depth + 1) * kr;
				break;
			}
			default:
			{
				// [comment]
				// We use the Phong illumation model int the default case. The phong model
				// is composed of a diffuse and a specular reflection component.
				// [/comment]
				Vec3f lightAmt = 0, specularColor = 0;
				Vec3f shadowPointOrig = (dir.dotProduct(N) < 0) ?
					hitPoint + N * options.bias :
					hitPoint - N * options.bias;
				// [comment]
				// Loop over all lights in the scene and sum their contribution up
				// We also apply the lambert cosine law here though we haven't explained yet what this means.
				// [/comment]
				for (uint32_t i = 0; i < lights.size(); ++i) {
					state.numShadowRays++;
					Vec3f lightDir = lights[i]->position - hitPoint;
					// square of the distance between hitPoint and the light
					float lightDistance2 = lightDir.dotProduct(lightDir);
					lightDir = lightDir.normalize();
					float LdotN = std::max(0.f, lightDir.dotProduct(N));
					Object *shadowHitObject = nullptr;
					float tNearShadow = kInfinity;
					// is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
					bool inShadow = trace(shadowPointOrig, lightDir, objects, tNearShadow, index, uv, &shadowHitObject) &&
						tNearShadow * tNearShadow < lightDistance2;
					lightAmt += (1.f - inShadow) * lights[i]->intensity * LdotN;
					Vec3f reflectionDirection = reflect(-lightDir, N);
					specularColor += powf(std::max(0.f, -reflectionDirection.dotProduct(dir)), hitObject->specularExponent) * lights[i]->intensity;
				}
				hitColor = lightAmt * hitObject->evalDiffuseColor(st) * hitObject->Kd + specularColor * hitObject->Ks;
				break;
			}
		}
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
	const std::vector<std::unique_ptr<LightLite>> &lights,
	State &state)
{
	Vec3f *framebuffer = new Vec3f[options.width * options.height];
	Vec3f *pix = framebuffer;
	float scale = tan(scratch::utils::deg2rad(options.fov * 0.5f));
	float imageAspectRatio = options.width / (float)options.height;
	Vec3f orig(0);
	for (uint32_t j = 0; j < options.height; ++j) {
		for (uint32_t i = 0; i < options.width; ++i) {
			// generate primary ray direction
			float x = (2 * (i + 0.5f) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5f) / (float)options.height) * scale;
			Vec3f dir = Vec3f(x, y, -1).normalize();
			state.numPrimaryRays++;
			*(pix++) = castRay(orig, dir, objects, lights, options, state, 0);
		}
	}

	// save framebuffer to file
	std::ofstream ofs;
	ofs.open("./ray_whitted.ppm", std::ios::out | std::ios::binary);
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
	std::vector<std::unique_ptr<LightLite>> lights;

	Sphere *sph1 = new Sphere(Vec3f(-1.f, 0.f, -12.f), 2.f);
	sph1->materialType = kDiffuseAndGlossy;
	sph1->diffuseColor = Vec3f(0.6f, 0.7f, 0.8f);
	Sphere *sph2 = new Sphere(Vec3f(0.5f, -0.5f, -8.f), 1.5f);
	sph2->ior = 1.5;
	sph2->materialType = kReflectionAndRefraction;

	objects.push_back(std::unique_ptr<Sphere>(sph1));
	objects.push_back(std::unique_ptr<Sphere>(sph2));

	Vec3f verts[4] = {{-5.f,-3.f,-6.f},
					  {5.f,-3.f,-6.f},
					  {5.f,-3.f,-16.f},
					  {-5.f,-3.f,-16.f}};
	uint32_t vertIndex[6] = {0, 1, 3, 1, 2, 3};
	Vec2f st[4] = {{0.f, 0.f}, {1.f, 0.f}, {1.f, 1.f}, {0.f, 1.f}};
	TriangleMesh *mesh = new TriangleMesh(verts, vertIndex, 2, st);
	mesh->materialType = kDiffuseAndGlossy;

	objects.push_back(std::unique_ptr<TriangleMesh>(mesh));

	lights.push_back(std::unique_ptr<LightLite>(new LightLite(Vec3f(-20.f, 70.f, 20.f), 0.5f)));
	lights.push_back(std::unique_ptr<LightLite>(new LightLite(Vec3f(30.f, 50.f, -12.f), 1.f)));

	// setting up options
	Options options;
	options.width = 1920;
	options.height = 1080;
	options.fov = 90.f;
	options.backgroundColor = Vec3f(0.235294f, 0.67451f, 0.843137f);
	options.maxDepth = 50;
	options.bias = 0.00001f;

	//setting up state
	State state;
	state.numPrimaryRays = 0;
	state.numReflectionRays = 0;
	state.numRefractionRays = 0;
	state.numShadowRays = 0;

	auto t_start = std::chrono::high_resolution_clock::now();

	// finally, render
	render(options, objects, lights, state);

	auto t_end = std::chrono::high_resolution_clock::now();
	auto passedTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
	std::cerr << "Wall passed time:  " << passedTime << " ms" << std::endl;

	std::cerr << "Primary rays cast: " << state.numPrimaryRays << std::endl;
	std::cerr << "Reflection rays cast: " << state.numReflectionRays << std::endl;
	std::cerr << "Refraction rays cast: " << state.numRefractionRays << std::endl;
	std::cerr << "Shadow rays cast: " << state.numShadowRays << std::endl;

	return 0;
}
