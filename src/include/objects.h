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

#pragma once

#include <memory>

#include "geometry.h"

enum MaterialType {
	DIFFUSE_AND_GLOSSY,
	REFLECTION_AND_REFRACTION,
	REFLECTION
};

class Object
{
 public:
	Object() :
		materialType(DIFFUSE_AND_GLOSSY),
		ior(1.3f), Kd(0.8f), Ks(0.2f), diffuseColor(0.2f), specularExponent(25) {}
	virtual ~Object() {}
	virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;
	virtual void getSurfaceProperties(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec3f &, Vec2f &) const = 0;
	virtual Vec3f evalDiffuseColor(const Vec2f &) const { return diffuseColor; }
	// material properties
	MaterialType materialType;
	float ior;
	float Kd, Ks;
	Vec3f diffuseColor;
	float specularExponent;
};

class Sphere : public Object
{
public:
	Sphere(const Vec3f &c,
		   const float &r,
		   const Vec3f &sc = 0,
		   const float &refl = 0,
		   const float &transp = 0,
		   const Vec3f &ec = 0) :
		   center(c),
		   radius(r),
		   radius2(r * r),
		   surfaceColor(sc),
		   emissionColor(ec),
		   transparency(transp),
		   reflection(refl) {}

	//[comment]
	// Compute a ray-sphere intersection using the analytic solution
	//[/comment]
	bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, uint32_t &index, Vec2f &uv) const
	{
		// analytic solution
		Vec3f L = orig - center;
		float a = dir.dotProduct(dir);
		float b = 2 * dir.dotProduct(L);
		float c = L.dotProduct(L) - radius2;
		float t0, t1;
		if (!solveQuadratic(a, b, c, t0, t1)) return false;
		if (t0 < 0) t0 = t1;
		if (t0 < 0) return false;
		tnear = t0;

		return true;
	}

	/*
	//[comment]
	// Compute a ray-sphere intersection using the geometric solution
	//[/comment]
	bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
	{
			Vec3f l = center - rayorig;
			float tca = l.dotProduct(raydir);
			if (tca < 0) return false;
			float d2 = l.dotProduct(l) - tca * tca;
			if (d2 > radius2) return false;
			float thc = sqrt(radius2 - d2);
			t0 = tca - thc;
			t1 = tca + thc;

			return true;
	}
	*/

	void getSurfaceProperties(const Vec3f &P, const Vec3f &I, const uint32_t &index, const Vec2f &uv, Vec3f &N, Vec2f &st) const
	{
		N = (P - center).normalize();
	}

	Vec3f center;                           /// position of the sphere
	float radius, radius2;                  /// sphere radius and radius^2
	Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
	float transparency, reflection;         /// surface transparency and reflectivity
};

bool rayTriangleIntersect(
	const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
	const Vec3f &orig, const Vec3f &dir,
	float &tnear, float &u, float &v)
{
	Vec3f edge1 = v1 - v0;
	Vec3f edge2 = v2 - v0;
	Vec3f pvec = dir.crossProduct(edge2);
	float det = edge1.dotProduct(pvec);
	if (det == 0 || det < 0) return false;

	Vec3f tvec = orig - v0;
	u = tvec.dotProduct(pvec);
	if (u < 0 || u > det) return false;

	Vec3f qvec = tvec.crossProduct(edge1);
	v = dir.dotProduct(qvec);
	if (v < 0 || u + v > det) return false;

	float invDet = 1 / det;

	tnear = edge2.dotProduct(qvec) * invDet;
	u *= invDet;
	v *= invDet;

	return true;
}

class MeshTriangle : public Object
{
public:
	MeshTriangle(
		const Vec3f *verts,
		const uint32_t *vertsIndex,
		const uint32_t &numTris,
		const Vec2f *st)
	{
		uint32_t maxIndex = 0;
		for (uint32_t i = 0; i < numTris * 3; ++i)
			if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
		maxIndex += 1;
		vertices = std::unique_ptr<Vec3f[]>(new Vec3f[maxIndex]);
		memcpy(vertices.get(), verts, sizeof(Vec3f) * maxIndex);
		vertexIndex = std::unique_ptr<uint32_t[]>(new uint32_t[numTris * 3]);
		memcpy(vertexIndex.get(), vertsIndex, sizeof(uint32_t) * numTris * 3);
		numTriangles = numTris;
		stCoordinates = std::unique_ptr<Vec2f[]>(new Vec2f[maxIndex]);
		memcpy(stCoordinates.get(), st, sizeof(Vec2f) * maxIndex);
	}

	bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, uint32_t &index, Vec2f &uv) const
	{
		bool intersect = false;
		for (uint32_t k = 0; k < numTriangles; ++k) {
			const Vec3f & v0 = vertices[vertexIndex[k * 3]];
			const Vec3f & v1 = vertices[vertexIndex[k * 3 + 1]];
			const Vec3f & v2 = vertices[vertexIndex[k * 3 + 2]];
			float t, u, v;
			if (rayTriangleIntersect(v0, v1, v2, orig, dir, t, u, v) && t < tnear) {
				tnear = t;
				uv.x = u;
				uv.y = v;
				index = k;
				intersect |= true;
			}
		}

		return intersect;
	}

	void getSurfaceProperties(const Vec3f &P, const Vec3f &I, const uint32_t &index, const Vec2f &uv, Vec3f &N, Vec2f &st) const
	{
		const Vec3f &v0 = vertices[vertexIndex[index * 3]];
		const Vec3f &v1 = vertices[vertexIndex[index * 3 + 1]];
		const Vec3f &v2 = vertices[vertexIndex[index * 3 + 2]];
		Vec3f e0 = (v1 - v0).normalize();
		Vec3f e1 = (v2 - v1).normalize();
		N = e0.crossProduct(e1).normalize();
		const Vec2f &st0 = stCoordinates[vertexIndex[index * 3]];
		const Vec2f &st1 = stCoordinates[vertexIndex[index * 3 + 1]];
		const Vec2f &st2 = stCoordinates[vertexIndex[index * 3 + 2]];
		st = st0 * (1 - uv.x - uv.y) + st1 * uv.x + st2 * uv.y;
	}

	Vec3f evalDiffuseColor(const Vec2f &st) const
	{
		float scale = 5.f;
		float pattern = (fmodf(st.x * scale, 1.f) > 0.5f) ^ (fmodf(st.y * scale, 1.f) > 0.5f);
		return mix(Vec3f(0.815f, 0.235f, 0.031f), Vec3f(0.937f, 0.937f, 0.231f), pattern);
	}

	std::unique_ptr<Vec3f[]> vertices;
	uint32_t numTriangles;
	std::unique_ptr<uint32_t[]> vertexIndex;
	std::unique_ptr<Vec2f[]> stCoordinates;
};

class Light
{
public:
	Light(const Vec3f &p, const Vec3f &i) : position(p), intensity(i) {}
	Vec3f position;
	Vec3f intensity;
};