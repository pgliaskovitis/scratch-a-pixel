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
// Example of Acceleration Structures for Ray-Tracing (BBox, BVH and Grid)
//[/header]

#include <iostream>
#include <fstream>
#include <chrono>

#include "objects.h"
#include "acceleration.h"
#include "assert.h"
#include "teapot_data.h"

struct Options
{
	float fov;
	uint32_t width;
	uint32_t height;
	Matrix44f cameraToWorld;
};

std::vector<std::unique_ptr<const TriangleMesh>> createUtahTeapot(const Options& options)
{
	Matrix44f rotate90(1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1);
	std::vector<std::unique_ptr<const TriangleMesh>> meshes;
	uint32_t width = 8, height = 8;
	uint32_t numPolygons = width * height;
	std::unique_ptr<uint32_t []> polyNumVertsArray(new uint32_t[numPolygons]);
	for (uint32_t i = 0; i < numPolygons; ++i) {
		polyNumVertsArray[i] = 4;
	}
	std::unique_ptr<uint32_t []> polyIndicesInVertPool(new uint32_t[numPolygons * 4]);
	// set indices
	for (uint32_t y = 0, offset = 0; y < height; ++y) {
		for (uint32_t x = 0; x < width; ++x, offset += 4) {
			// counter-clockwise to get the normal pointing in the right direction
			polyIndicesInVertPool[offset    ] = (width + 1) * y + x;
			polyIndicesInVertPool[offset + 3] = (width + 1) * y + x + 1;
			polyIndicesInVertPool[offset + 2] = (width + 1) * (y + 1) + x + 1;
			polyIndicesInVertPool[offset + 1] = (width + 1) * (y + 1) + x;
		}
	}
	Vec3f controlPoints[16];
	for (uint32_t i  = 0; i < kTeapotNumPatches; ++i) {
		std::unique_ptr<Vec3f []> vertPool(new Vec3f[(width + 1) * (height + 1)]);
		std::unique_ptr<Vec3f []> N(new Vec3f[(width + 1) * (height + 1)]);
		std::unique_ptr<Vec2f []> st(new Vec2f[(width + 1) * (height + 1)]);
		for (uint32_t j = 0; j < 16; ++j) {
			controlPoints[j].x = teapotVertices[teapotPatches[i][j] - 1][0],
			controlPoints[j].y = teapotVertices[teapotPatches[i][j] - 1][1],
			controlPoints[j].z = teapotVertices[teapotPatches[i][j] - 1][2];
		}
		for (uint32_t y = 0, currVertIndex = 0; y <= height; ++y) {
			float v = y / (float)height;
			for (uint32_t x = 0; x <= width; ++x, ++currVertIndex) {
				float u = x / (float)width;
				vertPool[currVertIndex] = evalBezierPatch(controlPoints, u, v);
				rotate90.multDirMatrix(vertPool[currVertIndex], vertPool[currVertIndex]);
				Vec3f dU = dUBezier(controlPoints, u, v);
				Vec3f dV = dVBezier(controlPoints, u, v);
				N[currVertIndex] = dU.crossProduct(dV).normalize();
				st[currVertIndex].x = u;
				st[currVertIndex].y = v;
			}
		}

		meshes.emplace_back(new TriangleMesh(Matrix44f::kIdentity, numPolygons, polyNumVertsArray, polyIndicesInVertPool, vertPool, N, st));
	}

	return meshes;
}

void makeScene(std::vector<std::unique_ptr<const TriangleMesh>>& meshes, const Options& options)
{
	meshes = std::move(createUtahTeapot(options));
}

// [comment]
// Main Render() function. Loop over each pixel in the image and trace a primary
// ray starting from the camera origin and passing through the current pixel. If they
// ray intersects geometry in the scene return some color (the color of the object)
// otherwise nothing (black or the color of the background)
// [/comment]
void render(const std::unique_ptr<AccelerationStructure>& accel, const Options& options)
{
	std::unique_ptr<Vec3f []> buffer(new Vec3f[options.width * options.height]);
	Vec3f orig(0.f, 0.f, 5.f);
	options.cameraToWorld.multVecMatrix(orig, orig);
	float scale = std::tan(scratch::utils::deg2rad(options.fov * 0.5));
	float imageAspectRatio = options.width / static_cast<float>(options.height);
	assert(imageAspectRatio > 1);
	uint32_t rayId = 1; // Start at 1 not 0!! (see Grid code and mailboxing)
	auto timeStart = std::chrono::high_resolution_clock::now();
	for (uint32_t j = 0; j < options.height; ++j) {
		for (uint32_t i = 0; i < options.width; ++i) {
			Vec3f dir((2 * (i + 0.5f) / options.width - 1) * scale * imageAspectRatio,
					  (1 - 2 * (j + 0.5) / options.height) * scale, -1);
			options.cameraToWorld.multDirMatrix(dir, dir);
			dir.normalize();
			accel->getStats().numPrimaryRays++;
			float tHit = kInfinity;
			buffer[j * options.width + i] = (accel->intersect(orig, dir, rayId++, tHit)) ? Vec3f(1) :  Vec3f(0);
			// fprintf(stderr, "\r%3d%c", uint32_t(j / (float)options.height * 100), '%');
		}
	}
	auto timeEnd = std::chrono::high_resolution_clock::now();
	auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
	fprintf(stderr, "\rDone: %.2f (sec)\n", passedTime / 1000);

	// store to PPM file
	std::ofstream ofs;
	ofs.open("ray_acceleration.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.width * options.height; ++i) {
		Vec3<uint8_t> pixRgb;
		pixRgb.x = static_cast<uint8_t>(255 * std::max(0.f, std::min(1.f, buffer[i].x)));
		pixRgb.y = static_cast<uint8_t>(255 * std::max(0.f, std::min(1.f, buffer[i].y)));
		pixRgb.z = static_cast<uint8_t>(255 * std::max(0.f, std::min(1.f, buffer[i].z)));
		ofs << pixRgb.x << pixRgb.y << pixRgb.z;
	}
	ofs.close();
}

void exportTriangleMesh(const std::vector<std::unique_ptr<const TriangleMesh>>& meshes)
{
	std::ofstream f;
	f.open("mesh.obj");
	uint32_t k = 0, off = 0;
	for (const auto& mesh : meshes) {
		f << "g default" << std::endl;
		for (uint32_t i = 0; i < mesh->numTris * 3; ++i) {
			f << "v " << mesh->P[i].x << " " << mesh->P[i].y << " " << mesh->P[i].z << std::endl;
		}

		f << "g mesh" << k++ << std::endl;
		for (uint32_t i = 0; i < mesh->numTris; ++i) {
			f << "f " << mesh->trisIndex[i * 3] + 1 + off << " " <<
			mesh->trisIndex[i * 3 + 1] + 1 + off << " " <<
			mesh->trisIndex[i * 3 + 2] + 1 + off << std::endl;
		}
		off += mesh->numTris * 3;
	}

	f.close();
}

int main(int argc, char **argv)
{
	Options options;
	options.width = 1920;
	options.height = 1080;
	options.fov = 90.0f;

	std::vector<std::unique_ptr<const TriangleMesh>> meshes;
	makeScene(meshes, options);
	// exportTriangleMesh(meshes);
	uint32_t numTriangles = 0;
	for (const auto& mesh : meshes) {
		numTriangles += mesh->numTris;
	}

	// [comment]
	// Create the acceleration structure
	// [/comment]
// #if defined(ACCEL_BBOX)
//	std::unique_ptr<AccelerationStructure> accel(new BBoxAcceleration(meshes));
// #elif defined(ACCEL_BVH)
//	std::unique_ptr<AccelerationStructure> accel(new BVH(meshes));
// #elif defined(ACCEL_GRID)
	std::unique_ptr<AccelerationStructure> accel(new Grid(meshes));
// #else
//	std::unique_ptr<AccelerationStructure> accel(new AccelerationStructure(meshes));
// #endif

	render(accel, options);

	std::cout << "Total number of triangles                         | " << numTriangles << std::endl;
	std::cout << "Total number of primary rays                      | " << accel->getStats().numPrimaryRays << std::endl;
	std::cout << "Total number of ray-mesh tests                    | " << accel->getStats().numRayMeshTests << std::endl;
	std::cout << "Total number of ray-bbox tests                    | " << accel->getStats().numRayBBoxTests << std::endl;
	std::cout << "Total number of ray-bounding volume tests         | " << accel->getStats().numRayBoundingVolumeTests << std::endl;
	std::cout << "Total number of ray-bounding volume intersections | " << accel->getStats().numRayBoundingVolumeIntersections << std::endl;
	std::cout << "Total number of ray-triangle tests                | " << accel->getStats().numRayTriangleTests << std::endl;
	std::cout << "Total number of ray-triangle intersections        | " << accel->getStats().numRayTriangleIntersections << std::endl;

	return 0;
}
