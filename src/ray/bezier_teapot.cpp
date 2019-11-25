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
// This program generate and render the Utah teapot
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
#include "utils.h"
#include "lights.h"
#include "loader.h"

#include "teapot_data.h"
#include "hair_data.h"

struct Options
{
	uint32_t width = 1920;
	uint32_t height = 1080;
	float fov = 90;
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
	IntersectInfo &isect)
{
	isect.hitObject = nullptr;
	for (uint32_t k = 0; k < objects.size(); ++k) {
		float tNearTriangle = kInfinity;
		uint32_t indexTriangle;
		Vec2f uvTriangle;
		if (objects[k]->intersect(orig, dir, tNearTriangle, indexTriangle, uvTriangle) && tNearTriangle < isect.tNear) {
			isect.hitObject = objects[k].get();
			isect.tNear = tNearTriangle;
			isect.index = indexTriangle;
			isect.uv = uvTriangle;
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
	if (depth > options.maxDepth) return 0;
	Vec3f hitColor = 0;
	IntersectInfo isect;
	if (trace(orig, dir, objects, isect)) {
		Vec3f hitPoint = orig + dir * isect.tNear;
		Vec3f hitNormal;
		Vec2f hitTexCoordinates;
		isect.hitObject->getSurfaceProperties(hitPoint, dir, isect.index, isect.uv, hitNormal, hitTexCoordinates);
		hitColor = std::max(0.f, -hitNormal.dotProduct(dir)) ;//* Vec3f(hitTexCoordinates.x, hitTexCoordinates.y, 1);
	}
	else {
		hitColor = 0.3f;
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
	Vec3f *framebuffer  = new Vec3f[options.width * options.height];
	Vec3f *pix = framebuffer;
	float scale = tan(scratch::utils::deg2rad(options.fov * 0.5));
	float imageAspectRatio = options.width / (float)options.height;
	Vec3f orig;
	options.cameraToWorld.multVecMatrix(Vec3f(0), orig);
	auto timeStart = std::chrono::high_resolution_clock::now();
	for (uint32_t j = 0; j < options.height; ++j) {
		for (uint32_t i = 0; i < options.width; ++i) {
			// generate primary ray direction
			float x = (2 * (i + 0.5) / (float)options.width - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5) / (float)options.height) * scale;
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
	ofs.open("ray_bezier_teapot.ppm");
	ofs << "P6\n" << options.width << " " << options.height << "\n255\n";
	for (uint32_t i = 0; i < options.height * options.width; ++i) {
		char r = (char)(255 * scratch::utils::clamp(0, 1, powf(framebuffer[i].x, 1 / gamma)));
		char g = (char)(255 * scratch::utils::clamp(0, 1, powf(framebuffer[i].y, 1 / gamma)));
		char b = (char)(255 * scratch::utils::clamp(0, 1, powf(framebuffer[i].z, 1 / gamma)));
		ofs << r << g << b;
	}
	ofs.close();
	delete [] framebuffer;
}

// [comment]
// Generate a poly-mesh Utah teapot out of Bezier patches
// [/comment]
void createPolyTeapot(const Matrix44f& o2w, std::vector<std::unique_ptr<Object>> &objects)
{
	uint32_t divs = 8;
	std::unique_ptr<Vec3f []> P(new Vec3f[(divs + 1) * (divs + 1)]);
	std::unique_ptr<uint32_t []> nvertices(new uint32_t[divs * divs]);
	std::unique_ptr<uint32_t []> vertices(new uint32_t[divs * divs * 4]);
	std::unique_ptr<Vec3f []> N(new Vec3f[(divs + 1) * (divs + 1)]);
	std::unique_ptr<Vec2f []> st(new Vec2f[(divs + 1) * (divs + 1)]);

	// face connectivity - all patches are subdivided the same way so there
	// share the same topology and uvs
	for (uint16_t j = 0, k = 0; j < divs; ++j) {
		for (uint16_t i = 0; i < divs; ++i, ++k) {
			nvertices[k] = 4;
			vertices[k * 4    ] = (divs + 1) * j + i;
			vertices[k * 4 + 1] = (divs + 1) * j + i + 1;
			vertices[k * 4 + 2] = (divs + 1) * (j + 1) + i + 1;
			vertices[k * 4 + 3] = (divs + 1) * (j + 1) + i;
		}
	}

	Vec3f controlPoints[16];
	for (int np = 0; np < kTeapotNumPatches; ++np) { // kTeapotNumPatches
		// set the control points for the current patch
		for (uint32_t i = 0; i < 16; ++i) {
			controlPoints[i][0] = teapotVertices[teapotPatches[np][i] - 1][0],
			controlPoints[i][1] = teapotVertices[teapotPatches[np][i] - 1][1],
			controlPoints[i][2] = teapotVertices[teapotPatches[np][i] - 1][2];

			// generate grid
			for (uint16_t j = 0, k = 0; j <= divs; ++j) {
				float v = j / (float)divs;
				for (uint16_t i = 0; i <= divs; ++i, ++k) {
					float u = i / (float)divs;
					P[k] = evalBezierPatch(controlPoints, u, v);
					Vec3f dU = dUBezier(controlPoints, u, v);
					Vec3f dV = dVBezier(controlPoints, u, v);
					N[k] = dU.crossProduct(dV).normalize();
					st[k].x = u;
					st[k].y = v;
				}
			}
		}

		auto triangle_mesh = std::make_unique<TriangleMesh>(TriangleMesh(o2w, divs * divs, nvertices, vertices, P, N, st));
		triangle_mesh->specularExponent = 10;
		objects.push_back(std::move(triangle_mesh));
	}
}

// [comment]
// Generate a thin cylinder centred around a Bezier curve
// [/comment]
void createCurveGeometry(std::vector<std::unique_ptr<Object>> &objects)
{
	uint32_t ndivs = 16;
	uint32_t ncurves = 1 + (curveNumPts - 4) / 3;
	Vec3f pts[4];
	std::unique_ptr<Vec3f []> P(new Vec3f[(ndivs + 1) * ndivs * ncurves + 1]);
	std::unique_ptr<Vec3f []> N(new Vec3f[(ndivs + 1) * ndivs * ncurves + 1]);
	std::unique_ptr<Vec2f []> st(new Vec2f[(ndivs + 1) * ndivs * ncurves + 1]);
	for (uint32_t i = 0; i < ncurves; ++i) {
		for (uint32_t j = 0; j < ndivs; ++j) {
			pts[0] = curveData[i * 3];
			pts[1] = curveData[i * 3 + 1];
			pts[2] = curveData[i * 3 + 2];
			pts[3] = curveData[i * 3 + 3];
			float s = j / (float)ndivs;
			Vec3f pt = evalBezierCurve(pts, s);
			Vec3f tangent = derivBezier(pts, s).normalize();
			bool swap = false;

			uint8_t maxAxis;
			if (std::abs(tangent.x) > std::abs(tangent.y)) {
				if (std::abs(tangent.x) > std::abs(tangent.z)) {
					maxAxis = 0;
				}
				else {
					maxAxis = 2;
				}
			} else if (std::abs(tangent.y) > std::abs(tangent.z)) {
				maxAxis = 1;
			} else {
				maxAxis = 2;
			}
			Vec3f up, forward, right;

			switch (maxAxis) {
				case 0:
				case 1:
					up = tangent;
					forward = Vec3f(0, 0, 1);
					right = up.crossProduct(forward);
					forward = right.crossProduct(up);
					break;
				case 2:
					up = tangent;
					right = Vec3f(0, 0, 1);
					forward = right.crossProduct(up);
					right = up.crossProduct(forward);
					break;
				default:
					break;
			};

			float sNormalized = (i * ndivs + j) / float(ndivs * ncurves);
			float rad = 0.1 * (1 - sNormalized);
			for (uint32_t k = 0; k <= ndivs; ++k) {
				float t = k / (float)ndivs;
				float theta = t * 2 * M_PI;
				Vec3f pc(cos(theta) * rad, 0, sin(theta) * rad);
				float x = pc.x * right.x + pc.y * up.x + pc.z * forward.x;
				float y = pc.x * right.y + pc.y * up.y + pc.z * forward.y;
				float z = pc.x * right.z + pc.y * up.z + pc.z * forward.z;
				P[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vec3f(pt.x + x, pt.y + y, pt.z + z);
				N[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vec3f(x, y, z).normalize();
				st[i * (ndivs + 1) * ndivs + j * (ndivs + 1) + k] = Vec2f(sNormalized, t);
			}
		}
	}

	P[(ndivs + 1) * ndivs * ncurves] = curveData[curveNumPts - 1];
	N[(ndivs + 1) * ndivs * ncurves] = (curveData[curveNumPts - 2] - curveData[curveNumPts - 1]).normalize();
	st[(ndivs + 1) * ndivs * ncurves] = Vec2f(1, 0.5);
	uint32_t numFaces = ndivs * ndivs * ncurves;
	std::unique_ptr<uint32_t []> verts(new uint32_t[numFaces]);
	for (uint32_t i = 0; i < numFaces; ++i) {
		verts[i] = (i < (numFaces - ndivs)) ? 4 : 3;
	}
	std::unique_ptr<uint32_t []> vertIndices(new uint32_t[ndivs * ndivs * ncurves * 4 + ndivs * 3]);
	uint32_t nf = 0, ix = 0;
	for (uint32_t k = 0; k < ncurves; ++k) {
		for (uint32_t j = 0; j < ndivs; ++j) {
			if (k == (ncurves - 1) && j == (ndivs - 1)) { break; }
			for (uint32_t i = 0; i < ndivs; ++i) {
				vertIndices[ix] = nf;
				vertIndices[ix + 1] = nf + (ndivs + 1);
				vertIndices[ix + 2] = nf + (ndivs + 1) + 1;
				vertIndices[ix + 3] = nf + 1;
				ix += 4;
				++nf;
			}
			nf++;
		}
	}

	for (uint32_t i = 0; i < ndivs; ++i) {
		vertIndices[ix] = nf;
		vertIndices[ix + 1] = (ndivs + 1) * ndivs * ncurves;
		vertIndices[ix + 2] = nf + 1;
		ix += 3;
		nf++;
	}

	auto triangle_mesh = std::make_unique<TriangleMesh>(Matrix44f::kIdentity, numFaces, verts, vertIndices, P, N, st);
	triangle_mesh->specularExponent = 10;
	objects.push_back(std::move(triangle_mesh));
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

	createPolyTeapot(Matrix44f(1.f, 0.f, 0.f, 0.f, 0.f, 0.f, -1.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 1.f), objects);
	// createCurveGeometry(objects);

	// lights
	std::vector<std::unique_ptr<Light>> lights;
	Options options;

	// aliasing example
	options.fov = 39.89f;
	options.width = 1920;
	options.height = 1080;
	options.maxDepth = 1;

	// to render the teapot
	options.cameraToWorld = Matrix44f(0.897258f, 0.f, -0.441506f, 0.f, -0.288129f, 0.757698f, -0.585556f, 0.f, 0.334528f, 0.652606f, 0.679851f, 0.f, 5.439442f, 11.080794f, 10.381341f, 1.f);

	// to render the curve as geometry
	// options.cameraToWorld = Matrix44f(0.707107f, 0.f, -0.707107f, 0.f, -0.369866f, 0.85229f, -0.369866f, 0.f, 0.60266f, 0.523069f, 0.60266f, 0.f, 2.634f, 3.178036f, 2.262122f, 1.f);

	// finally, render
	render(options, objects, lights);

	return 0;
}
