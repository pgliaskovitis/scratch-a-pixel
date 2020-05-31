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

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <chrono>

#include "geometry.h"
#include "lights.h"
#include "implicit.h"

std::vector<std::shared_ptr<ImplicitShape>> makeScene()
{
	std::vector<std::shared_ptr<ImplicitShape>> shapes;

#if 0
	shapes.push_back(std::make_shared<ImplicitPlane>(Vec3f(0, 1, 0), Vec3f(0, -2, 0)));
	shapes.push_back(std::make_shared<ImplicitCube>(Vec3f(1.5)));
	shapes.push_back(std::make_shared<ImplicitTorus>(2, 0.65));
#elif 0
	shapes.push_back(std::make_shared<ImplicitPlane>(Vec3f(0, 1, 0), Vec3f(0, -2, 0)));
	shapes.push_back(std::make_shared<Blend>(
	std::make_shared<ImplicitCube>(Vec3f(1.5)),
	std::make_shared<ImplicitTorus>(2, 0.65), 5));
#elif 0
	shapes.push_back(std::make_shared<Blend>(
	std::make_shared<ImplicitPlane>(Vec3f(0, 1, 0), Vec3f(0, 0, 0)),
	std::make_shared<ImplicitTorus>(2, 0.65), 5));
#elif 0
	shapes.push_back(std::make_shared<ImplicitPlane>(Vec3f(0, 1, 0), Vec3f(0, -2, 0)));
	shapes.push_back(std::make_shared<Mix>(
	std::make_shared<ImplicitCube>(Vec3f(1)),
	std::make_shared<ImplicitSphere>(Vec3f(0), 1), 0.5));
#else
	shapes.push_back(std::make_shared<ImplicitPlane>(Vec3f(0, 1, 0), Vec3f(0, -2, 0)));
	shapes.push_back(std::make_shared<SoftObject>());
#endif
	return shapes;
}

bool sphereTraceShadow(
	const Vec3f& rayOrigin,
	const Vec3f& rayDirection,
	const float& maxDistance,
	const std::vector<std::shared_ptr<ImplicitShape>> scene)
{
	constexpr float threshold = 10e-5;
	float t = 0;

	while (t < maxDistance) {
		float minDistance = kInfinity;
		Vec3f from = rayOrigin + t * rayDirection;
		for (auto shape : scene) {
			float d = shape->getDistance(from);
			if (d < minDistance) {
				minDistance = d;
			}
		// did we find an intersection?
		if (minDistance <= threshold * t) {
			return true;
		}
	}

	// no intersection, move along the ray by minDistance
	t += minDistance;
}

	return false;
}

Vec3f shade(
	const Vec3f& rayOrigin,
	const Vec3f& rayDirection,
	const float& t,
	const ImplicitShape *shape,
	const std::vector<std::shared_ptr<ImplicitShape>> scene,
	const std::vector<std::unique_ptr<PointLight>> &lights)
{
	constexpr float delta = 10e-5;
	Vec3f p = rayOrigin + t * rayDirection;
	Vec3f n = Vec3f(
	shape->getDistance(p + Vec3f(delta, 0, 0)) - shape->getDistance(p + Vec3f(-delta, 0, 0)),
	shape->getDistance(p + Vec3f(0, delta, 0)) - shape->getDistance(p + Vec3f(0, -delta, 0)),
	shape->getDistance(p + Vec3f(0, 0, delta)) - shape->getDistance(p + Vec3f(0, 0, -delta)));
	n.normalize();

	Vec3f R = 0;

	// loop over all lights in the scene and add their contribution to P's brightness
	for (const auto& light: lights) {
		Vec3f lightDir = light->pos - p;
		if (lightDir.dotProduct(n) > 0) {
			float dist2 = lightDir.norm();
			lightDir.normalize();
			bool shadow = 1 - sphereTraceShadow(p, lightDir, sqrtf(dist2), scene);
			R += shadow * lightDir.dotProduct(n) * light->color * light->intensity / (4 * M_PI * dist2);
		}
	}

	return R;
}

Vec3f sphereTrace(
    const Vec3f& rayOrigin,
    const Vec3f& rayDirection,
    const std::vector<std::shared_ptr<ImplicitShape>>& scene,
    const std::vector<std::unique_ptr<PointLight>>& lights)
{
	constexpr float maxDistance = 100;
	float t = 0;
	uint32_t numSteps = 0;
	const ImplicitShape *isectShape = nullptr;

	constexpr float threshold = 10e-6;

	while (t < maxDistance) {
		float minDistance = kInfinity;
		Vec3f from = rayOrigin + t * rayDirection;
		for (const auto& shape : scene) {
			float d = shape->getDistance(from);
			if (d < minDistance) {
				minDistance = d;
				isectShape = shape.get();
			}
		}

		if (minDistance <= threshold * t) {
			return shade(rayOrigin, rayDirection, t, isectShape, scene, lights);
		}
		t += minDistance;
		numSteps++;
	}

	return 0;
}

int main(int argc, char **argv)
{
	srand48(13);

	Matrix44f camToWorld = lookAt(Vec3f(0, 1, 9), 0);

	std::vector<std::shared_ptr<ImplicitShape>> scene = makeScene();
	std::vector<std::unique_ptr<PointLight>> lights;
	lights.push_back(std::make_unique<PointLight>(Vec3f( 20, 30,  20), Vec3f(1.0, 0.9, 0.7), 4000));
	lights.push_back(std::make_unique<PointLight>(Vec3f(-20, 30, -20), Vec3f(0.8, 0.9, 1.0), 4000));
	lights.push_back(std::make_unique<PointLight>(Vec3f( -5, 10,  20), Vec3f(1.0, 1.0, 1.0), 3000));

	constexpr uint32_t width = 1920, height = 1080;
	constexpr float ratio = width / static_cast<float>(height);
	constexpr float fov = 60;
	float angle = tan(scratch::utils::deg2rad(fov * 0.5));
	Vec3f *buffer = new Vec3f[width * height];
	Vec3f rayOrigin;
	camToWorld.multVecMatrix(Vec3f(0), rayOrigin);

	auto timeStart = std::chrono::high_resolution_clock::now();
	for (uint32_t j = 0; j < height; ++j) {
		for (uint32_t i = 0; i < width; ++i) {
			float x = (2 * i / static_cast<float>(width) - 1) * ratio * angle;
			float y = (1 - j / static_cast<float>(height) * 2) * angle;
			Vec3f rayDirection;
			camToWorld.multDirMatrix(Vec3f(x, y, -1).normalize(), rayDirection);
			ImplicitShape *tmp;
			Vec3f pixelColor = sphereTrace(rayOrigin, rayDirection, scene, lights);
			buffer[width * j + i] = pixelColor;
		}
	}
	auto timeEnd = std::chrono::high_resolution_clock::now();
	auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
	fprintf(stderr, "\rDone: %.2f (sec)\n", passedTime / 1000);

	std::ofstream ofs;
	ofs.open("./ray_distance_fields.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (uint32_t i = 0; i < width * height; ++i) {
		unsigned char r = static_cast<unsigned char>(std::min(1.0f, buffer[i][0]) * 255);
		unsigned char g = static_cast<unsigned char>(std::min(1.0f, buffer[i][1]) * 255);
		unsigned char b = static_cast<unsigned char>(std::min(1.0f, buffer[i][2]) * 255);
		ofs << r << g << b;
	}
	ofs.close();

	delete [] buffer;

	return 0;
}