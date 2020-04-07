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

// [header]
// A simple implementation of Perlin and Improved Perlin Noise
// [/header]

#include <cmath>
#include <cstdio>
#include <random>
#include <functional>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "geometry.h"
#include "utils.h"

// if you use windows uncomment these two lines
//#include "stdafx.h"
//#define _USE_MATH_DEFINES

#define ANALYTICAL_NORMALS 1

class PerlinNoise
{
public:
	PerlinNoise(const unsigned &seed = 2020)
	{
		std::mt19937 generator(seed);
		std::uniform_real_distribution<float> distribution;
		auto dice = std::bind(distribution, generator);
		for (unsigned i = 0; i < tableSize; ++i) {
#if 0
			// bad
			float gradientLen2;
			do {
				gradients[i] = Vec3f(2 * dice() - 1, 2 * dice() - 1, 2 * dice() - 1);
				gradientLen2 = gradients[i].length2();
			} while (gradientLen2 > 1);
			gradients[i].normalize();
#else
			// better
			float theta = acos(2 * dice() - 1);
			float phi = 2 * dice() * M_PI;

			float x = cos(phi) * sin(theta);
			float y = sin(phi) * sin(theta);
			float z = cos(theta);
			gradients[i] = Vec3f(x, y, z);
#endif
			permutationTable[i] = i;
		}

		std::uniform_int_distribution<unsigned> distributionInt;
		auto diceInt = std::bind(distributionInt, generator);
		// create permutation table
		for (unsigned i = 0; i < tableSize; ++i) {
			std::swap(permutationTable[i], permutationTable[diceInt() & tableSizeMask]);
		}
		// extend the permutation table in the index range [256:512]
		for (unsigned i = 0; i < tableSize; ++i) {
			permutationTable[tableSize + i] = permutationTable[i];
		}
	}
	virtual ~PerlinNoise() {}

	//[comment]
	// Improved Noise implementation (2002)
	// This version compute the derivative of the noise function as well
	//[/comment]
	float eval(const Vec3f &p, Vec3f& derivs) const
	{
		int xi0 = ((int)std::floor(p.x)) & tableSizeMask;
		int yi0 = ((int)std::floor(p.y)) & tableSizeMask;
		int zi0 = ((int)std::floor(p.z)) & tableSizeMask;

		int xi1 = (xi0 + 1) & tableSizeMask;
		int yi1 = (yi0 + 1) & tableSizeMask;
		int zi1 = (zi0 + 1) & tableSizeMask;

		float tx = p.x - ((int)std::floor(p.x));
		float ty = p.y - ((int)std::floor(p.y));
		float tz = p.z - ((int)std::floor(p.z));

		float u = scratch::utils::quintic(tx);
		float v = scratch::utils::quintic(ty);
		float w = scratch::utils::quintic(tz);

		// generate vectors going from the grid points to p
		float x0 = tx, x1 = tx - 1;
		float y0 = ty, y1 = ty - 1;
		float z0 = tz, z1 = tz - 1;

		float a = gradientDotV(hash(xi0, yi0, zi0), x0, y0, z0);
		float b = gradientDotV(hash(xi1, yi0, zi0), x1, y0, z0);
		float c = gradientDotV(hash(xi0, yi1, zi0), x0, y1, z0);
		float d = gradientDotV(hash(xi1, yi1, zi0), x1, y1, z0);
		float e = gradientDotV(hash(xi0, yi0, zi1), x0, y0, z1);
		float f = gradientDotV(hash(xi1, yi0, zi1), x1, y0, z1);
		float g = gradientDotV(hash(xi0, yi1, zi1), x0, y1, z1);
		float h = gradientDotV(hash(xi1, yi1, zi1), x1, y1, z1);

		float du = scratch::utils::quinticDeriv(tx);
		float dv = scratch::utils::quinticDeriv(ty);
		float dw = scratch::utils::quinticDeriv(tz);

		float k0 = a;
		float k1 = (b - a);
		float k2 = (c - a);
		float k3 = (e - a);
		float k4 = (a + d - b - c);
		float k5 = (a + f - b - e);
		float k6 = (a + g - c - e);
		float k7 = (b + c + e + h - a - d - f - g);

		derivs.x = du *(k1 + k4 * v + k5 * w + k7 * v * w);
		derivs.y = dv *(k2 + k4 * u + k6 * w + k7 * v * w);
		derivs.z = dw *(k3 + k5 * u + k6 * v + k7 * v * w);

		return k0 + k1 * u + k2 * v + k3 * w + k4 * u * v + k5 * u * w + k6 * v * w + k7 * u * v * w;
	}

	//[comment]
	// classic/original Perlin noise implementation (1985)
	//[/comment]
	float eval(const Vec3f &p) const
	{
		int xi0 = ((int)std::floor(p.x)) & tableSizeMask;
		int yi0 = ((int)std::floor(p.y)) & tableSizeMask;
		int zi0 = ((int)std::floor(p.z)) & tableSizeMask;

		int xi1 = (xi0 + 1) & tableSizeMask;
		int yi1 = (yi0 + 1) & tableSizeMask;
		int zi1 = (zi0 + 1) & tableSizeMask;

		float tx = p.x - ((int)std::floor(p.x));
		float ty = p.y - ((int)std::floor(p.y));
		float tz = p.z - ((int)std::floor(p.z));

		float u = scratch::utils::smoothstep(tx);
		float v = scratch::utils::smoothstep(ty);
		float w = scratch::utils::smoothstep(tz);

		// gradients at the corner of the cell
		const Vec3f &c000 = gradients[hash(xi0, yi0, zi0)];
		const Vec3f &c100 = gradients[hash(xi1, yi0, zi0)];
		const Vec3f &c010 = gradients[hash(xi0, yi1, zi0)];
		const Vec3f &c110 = gradients[hash(xi1, yi1, zi0)];

		const Vec3f &c001 = gradients[hash(xi0, yi0, zi1)];
		const Vec3f &c101 = gradients[hash(xi1, yi0, zi1)];
		const Vec3f &c011 = gradients[hash(xi0, yi1, zi1)];
		const Vec3f &c111 = gradients[hash(xi1, yi1, zi1)];

		// generate vectors going from the grid points to p
		float x0 = tx, x1 = tx - 1;
		float y0 = ty, y1 = ty - 1;
		float z0 = tz, z1 = tz - 1;

		Vec3f p000 = Vec3f(x0, y0, z0);
		Vec3f p100 = Vec3f(x1, y0, z0);
		Vec3f p010 = Vec3f(x0, y1, z0);
		Vec3f p110 = Vec3f(x1, y1, z0);

		Vec3f p001 = Vec3f(x0, y0, z1);
		Vec3f p101 = Vec3f(x1, y0, z1);
		Vec3f p011 = Vec3f(x0, y1, z1);
		Vec3f p111 = Vec3f(x1, y1, z1);

		// linear interpolation
		float a = scratch::utils::lerp(c000.dotProduct(p000), c100.dotProduct(p100), u);
		float b = scratch::utils::lerp(c010.dotProduct(p010), c110.dotProduct(p110), u);
		float c = scratch::utils::lerp(c001.dotProduct(p001), c101.dotProduct(p101), u);
		float d = scratch::utils::lerp(c011.dotProduct(p011), c111.dotProduct(p111), u);

		float e = scratch::utils::lerp(a, b, v);
		float f = scratch::utils::lerp(c, d, v);

		return scratch::utils::lerp(e, f, w); // g
	}

private:
	/* inline */
	uint8_t hash(const int &x, const int &y, const int &z) const
	{
		return permutationTable[permutationTable[permutationTable[x] + y] + z];
	}

	//[comment]
	// Compute dot product between vector from cell corners to P with predefined gradient directions
	//
	//    perm: a value between 0 and 255
	//
	//    float x, float y, float z: coordinates of vector from cell corner to shaded point
	//[/comment]
	float gradientDotV(
		uint8_t perm, // a value between 0 and 255
		float x, float y, float z) const
	{
		switch (perm & 15) {
		case  0: return  x + y; // (1,1,0)
		case  1: return -x + y; // (-1,1,0)
		case  2: return  x - y; // (1,-1,0)
		case  3: return -x - y; // (-1,-1,0)
		case  4: return  x + z; // (1,0,1)
		case  5: return -x + z; // (-1,0,1)
		case  6: return  x - z; // (1,0,-1)
		case  7: return -x - z; // (-1,0,-1)
		case  8: return  y + z; // (0,1,1),
		case  9: return -y + z; // (0,-1,1),
		case 10: return  y - z; // (0,1,-1),
		case 11: return -y - z; // (0,-1,-1)
		case 12: return  y + x; // (1,1,0)
		case 13: return -x + y; // (-1,1,0)
		case 14: return -y + z; // (0,-1,1)
		case 15: return -y - z; // (0,-1,-1)
		}
	}

	static const unsigned tableSize = 256;
	static const unsigned tableSizeMask = tableSize - 1;
	Vec3f gradients[tableSize];
	unsigned permutationTable[tableSize * 2];
};

//[comment]
// Simple class to define a polygonal mesh
//[/comment]
class PolyMesh
{
public:
    PolyMesh() : vertices(nullptr), st(nullptr), normals(nullptr) {}
    ~PolyMesh()
    {
        if (vertices) delete[] vertices;
        if (st) delete[] st;
        if (normals) delete[] normals;
    }
    Vec3f *vertices;
    Vec2f *st;
    Vec3f *normals;
    uint32_t *faceArray;
    uint32_t *verticesArray;
    uint32_t numVertices;
    uint32_t numFaces;
    void exportToObj();
};

//[comment]
// Export polygonal mesh to OBJ file (vertex positions, texture coordinates and vertex normals)
//[/comment]
void PolyMesh::exportToObj()
{
	std::ofstream ofs;
	ofs.open("./noise_polymesh.obj", std::ios_base::out);

	for (uint32_t i = 0; i < numVertices; ++i) {
		ofs << "v " << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << std::endl;
	}

	for (uint32_t i = 0; i < numVertices; ++i) {
		ofs << "vt " << st[i].x << " " << st[i].y << std::endl;
	}

	for (uint32_t i = 0; i < numVertices; ++i) {
		ofs << "vn " << normals[i].x << " " << normals[i].y << " " << normals[i].z << std::endl;
	}

	for (uint32_t i = 0, k = 0; i < numFaces; ++i) {
		ofs << "f ";
		for (uint32_t j = 0; j < faceArray[i]; ++j) {
			uint32_t objIndex = verticesArray[k + j] + 1;
			ofs << objIndex << "/" << objIndex << "/" << objIndex << ((j == (faceArray[i] - 1)) ? "" : " ");
		}
		ofs << std::endl;
		k += faceArray[i];
	}

	ofs.close();
}

//[comment]
// Simple function to create a polygonal grid centred around the origin
//[/comment]
PolyMesh* createPolyMesh(
	uint32_t width = 1,
	uint32_t height = 1,
	uint32_t subdivisionWidth = 40,
	uint32_t subdivisionHeight = 40)
{
	PolyMesh *poly = new PolyMesh;
	poly->numVertices = (subdivisionWidth + 1) * (subdivisionHeight + 1);
	std::cerr << poly->numVertices << " vertices" << std::endl;
	poly->vertices = new Vec3f[poly->numVertices];
	poly->normals = new Vec3f[poly->numVertices];
	poly->st = new Vec2f[poly->numVertices];
	float invSubdivisionWidth = 1.f / subdivisionWidth;
	float invSubdivisionHeight = 1.f / subdivisionHeight;
	for (uint32_t j = 0; j <= subdivisionHeight; ++j) {
		for (uint32_t i = 0; i <= subdivisionWidth; ++i) {
			poly->vertices[j * (subdivisionWidth + 1) + i] = Vec3f(width * (i * invSubdivisionWidth - 0.5), 0, height * (j * invSubdivisionHeight - 0.5));
			poly->st[j * (subdivisionWidth + 1) + i] = Vec2f(i * invSubdivisionWidth, j * invSubdivisionHeight);
		}
	}

	poly->numFaces = subdivisionWidth * subdivisionHeight;
	poly->faceArray = new uint32_t[poly->numFaces];
	for (uint32_t i = 0; i < poly->numFaces; ++i) {
		poly->faceArray[i] = 4;
	}

	poly->verticesArray = new uint32_t[4 * poly->numFaces];
	for (uint32_t j = 0, k = 0; j < subdivisionHeight; ++j) {
		for (uint32_t i = 0; i < subdivisionWidth; ++i) {
			poly->verticesArray[k] = j * (subdivisionWidth + 1) + i;
			poly->verticesArray[k + 1] = j * (subdivisionWidth + 1) + i + 1;
			poly->verticesArray[k + 2] = (j + 1) * (subdivisionWidth + 1) + i + 1;
			poly->verticesArray[k + 3] = (j + 1) * (subdivisionWidth + 1) + i;
			k += 4;
		}
	}

	return poly;
}

int main(int argc, char **argv)
{
	PerlinNoise noise;

	PolyMesh *poly = createPolyMesh(3, 3, 900, 900);

	// displace and compute analytical normal using noise function partial derivatives
	Vec3f derivs;
	for (uint32_t i = 0; i < poly->numVertices; ++i) {
		Vec3f p((poly->vertices[i].x + 0.5), 0, (poly->vertices[i].z + 0.5));
		poly->vertices[i].y = noise.eval(p, derivs);
#if ANALYTICAL_NORMALS
		Vec3f tangent(1, derivs.x, 0); // tangent
		Vec3f bitangent(0, derivs.z, 1); // bitangent
		// equivalent to bitangent.cross(tangent)
		poly->normals[i] = Vec3f(-derivs.x, 1, -derivs.z);
		poly->normals[i].normalize();
#endif
	}

#if !ANALYTICAL_NORMALS
	// compute face normal if you want
	for (uint32_t k = 0, off = 0; k < poly->numFaces; ++k) {
		uint32_t nverts = poly->faceArray[k];
		const Vec3f &va = poly->vertices[poly->verticesArray[off]];
		const Vec3f &vb = poly->vertices[poly->verticesArray[off + 1]];
		const Vec3f &vc = poly->vertices[poly->verticesArray[off + nverts - 1]];

		Vec3f tangent = vb - va;
		Vec3f bitangent = vc - va;

		poly->normals[poly->verticesArray[off]] = bitangent.crossProduct(tangent);
		poly->normals[poly->verticesArray[off]].normalize();

		off += nverts;
	}
#endif

	poly->exportToObj();
	delete poly;

	// output noise map to PPM
	const uint32_t width = 2048, height = 2048;
	float *noiseMap = new float[width * height];

	for (uint32_t j = 0; j < height; ++j) {
		for (uint32_t i = 0; i < width; ++i) {
			noiseMap[j * width + i] = (noise.eval(Vec3f(i, 0, j) * (1 / 64.), derivs) + 1) * 0.5;
		}
	}

	std::ofstream ofs;
	ofs.open("./noise_perlin.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (unsigned k = 0; k < width * height; ++k) {
		unsigned char n = static_cast<unsigned char>(noiseMap[k] * 255);
		ofs << n << n << n;
	}
	ofs.close();

	delete[] noiseMap;

	return 0;
}