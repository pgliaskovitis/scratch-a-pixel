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

#include <atomic>
#include <cstring>
#include <random>

#include "geometry_utils.h"

enum MaterialType {
	kDiffuse,
	kDiffuseAndGlossy,
	kReflectionAndRefraction,
	kReflection,
	kPhong
};

class Ray
{
public:
	Ray(const Vec3f &orig, const Vec3f &dir) : orig(orig), dir(dir)
	{
		invdir = 1 / dir;
		sign[0] = (invdir.x < 0);
		sign[1] = (invdir.y < 0);
		sign[2] = (invdir.z < 0);
	}
	Vec3f orig, dir; // ray orig and dir
	Vec3f invdir;
	Vec3b sign;
};

template<typename T = float>
class BBox
{
public:
	BBox() {}

	BBox(Vec3<T> min_, Vec3<T> max_)
	{
		bounds[0] = min_;
		bounds[1] = max_;
	}

	BBox& extendBy(const Vec3<T>& p)
	{
		if (p.x < bounds[0].x) bounds[0].x = p.x;
		if (p.y < bounds[0].y) bounds[0].y = p.y;
		if (p.z < bounds[0].z) bounds[0].z = p.z;
		if (p.x > bounds[1].x) bounds[1].x = p.x;
		if (p.y > bounds[1].y) bounds[1].y = p.y;
		if (p.z > bounds[1].z) bounds[1].z = p.z;

		return *this;
	}

	Vec3<T> centroid() const { return (bounds[0] + bounds[1]) * 0.5; }
	Vec3<T>& operator [] (bool i) { return bounds[i]; }
	const Vec3<T> operator [] (bool i) const { return bounds[i]; }
	bool intersect(const Vec3<T>&, const Vec3<T>&, const Vec3b&, float&) const;
	bool intersect(const Ray &r, float &t) const;
	Vec3<T> bounds[2] = { kInfinity, -kInfinity };
};

template<typename T>
bool BBox<T>::intersect(const Vec3<T>& orig, const Vec3<T>& invDir, const Vec3b& sign, float& tHit) const
{
	float tmin, tmax, tymin, tymax, tzmin, tzmax;

	tmin  = (bounds[sign[0]    ].x - orig.x) * invDir.x;
	tmax  = (bounds[1 - sign[0]].x - orig.x) * invDir.x;
	tymin = (bounds[sign[1]    ].y - orig.y) * invDir.y;
	tymax = (bounds[1 - sign[1]].y - orig.y) * invDir.y;

	if ((tmin > tymax) || (tymin > tmax)) {
		return false;
	}

	if (tymin > tmin) {
		tmin = tymin;
	}
	if (tymax < tmax) {
		tmax = tymax;
	}

	tzmin = (bounds[sign[2]    ].z - orig.z) * invDir.z;
	tzmax = (bounds[1 - sign[2]].z - orig.z) * invDir.z;

	if ((tmin > tzmax) || (tzmin > tmax)) {
		return false;
	}

	if (tzmin > tmin) {
		tmin = tzmin;
	}
	if (tzmax < tmax) {
		tmax = tzmax;
	}

	tHit = tmin;

	/*
	if (tHit < 0) {
		tHit = tmax;
		if (tHit < 0) {
			return false;
		}
	}
	*/

	return true;
}

template <typename T>
bool BBox<T>::intersect(const Ray &r, float &t) const
{
	return intersect(r.orig, r.invdir, r.sign, t);
}

class Object
{
public:
	Object()
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0.f, 1.f);
		diffuseColor = Vec3f(dis(gen), dis(gen), dis(gen));
	}

	Object(const Matrix44f &o2w) : objectToWorld(o2w), worldToObject(o2w.inverse()) {}

	virtual ~Object() {}
	virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;
	virtual void getSurfaceProperties(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec3f &, Vec2f &) const = 0;
	virtual Vec3f evalDiffuseColor(const Vec2f &) const { return diffuseColor; }

	// material properties
	MaterialType materialType = kDiffuseAndGlossy;
	Vec3f diffuseColor; // albedo
	float ior = 1.3;
	float Kd = 0.8; // phong model diffuse weight
	float Ks = 0.2; // phong model specular weight
	float specularExponent = 25; // phong specular exponent

	// transforms
	Matrix44f objectToWorld;
	Matrix44f worldToObject;

	// shading
	bool smoothShading = false;

	// bounding box
	BBox<float> bbox;
};

class Sphere : public Object
{
public:

	Sphere(const Matrix44f &o2w, const float &r) : Object(o2w), radius(r), radius2(r * r)
	{
		o2w.multVecMatrix(Vec3f(0), center);
	}

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
		   reflection(refl)
	{}

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

	//[comment]
	// Compute a ray-sphere intersection using the analytic solution
	//[/comment]
	bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, uint32_t &index, Vec2f &uv) const
	{
		// analytic solution
		float t0, t1;
		Vec3f L = orig - center;
		float a = dir.dotProduct(dir);
		float b = 2 * dir.dotProduct(L);
		float c = L.dotProduct(L) - radius2;
		if (!scratch::utils::solveQuadratic(a, b, c, t0, t1)) return false;

		if (t0 < 0) {
			t0 = t1;
			if (t0 < 0) {
				return false;
			}
		}
		tnear = t0;

		return true;
	}

	void getSurfaceProperties(
		const Vec3f &hitPoint,
		const Vec3f &viewDirection,
		const uint32_t &index,
		const Vec2f &uv,
		Vec3f &hitNormal,
		Vec2f &hitTextureCoordinates) const
	{
		hitNormal = (hitPoint - center).normalize();
		// In this particular case, the normal is simular to a point on a unit sphere
		// centered around the origin. We can thus use the normal coordinates to compute
		// the spherical coordinates of Phit.
		// atan2 returns a value in the range [-pi, pi] and we need to remap it to range [0, 1]
		// acosf returns a value in the range [0, pi] and we also need to remap it to the range [0, 1]
		hitTextureCoordinates.x = (1 + atan2(hitNormal.z, hitNormal.x) / M_PI) * 0.5f;
		hitTextureCoordinates.y = acosf(hitNormal.y) / M_PI;
	}

	Vec3f center;                           /// position of the sphere
	float radius, radius2;                  /// sphere radius and radius^2
	Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
	float transparency, reflection;         /// surface transparency and reflectivity
};

class TriangleMesh : public Object
{
public:
	// Build a triangle mesh from a face index array and a vertex index array
	TriangleMesh(
		const Vec3f *verts,
		const uint32_t *vertsIndex,
		const uint32_t &numTris,
		const Vec2f *st)
	{
		uint32_t maxIndex = 0;
		for (uint32_t i = 0; i < numTris * 3; ++i) {
			if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
		}
		maxIndex += 1;
		P = std::unique_ptr<Vec3f[]>(new Vec3f[maxIndex]);
		memcpy(P.get(), verts, sizeof(Vec3f) * maxIndex);
		for (uint32_t i = 0; i < maxIndex; ++i) {
			bbox.extendBy(P[i]);
		}
		trisIndex = std::unique_ptr<uint32_t[]>(new uint32_t[numTris * 3]);
		memcpy(trisIndex.get(), vertsIndex, sizeof(uint32_t) * numTris * 3);
		this->numTris = numTris;
		texCoordinates = std::unique_ptr<Vec2f[]>(new Vec2f[maxIndex]);
		memcpy(texCoordinates.get(), st, sizeof(Vec2f) * maxIndex);
		mailbox.resize(numTris, 0); // should all be initialized with 0
	}

	// Build a triangle mesh from a face index array and a vertex index array
	TriangleMesh(
		const uint32_t nfaces,
		const std::unique_ptr<uint32_t []> &faceIndex,
		const std::unique_ptr<uint32_t []> &vertsIndex,
		const std::unique_ptr<Vec3f []> &verts,
		std::unique_ptr<Vec3f []> &normals,
		std::unique_ptr<Vec2f []> &st) :
			numTris(0)
	{
		uint32_t k = 0, maxVertIndex = 0;
		// find out how many triangles we need to create for this mesh
		for (uint32_t i = 0; i < nfaces; ++i) {
			numTris += faceIndex[i] - 2;
			for (uint32_t j = 0; j < faceIndex[i]; ++j) {
				if (vertsIndex[k + j] > maxVertIndex)
					maxVertIndex = vertsIndex[k + j];
			}
			k += faceIndex[i];
		}
		maxVertIndex += 1;

		// allocate memory to store the position of the mesh vertices
		P = std::unique_ptr<Vec3f []>(new Vec3f[maxVertIndex]);
		for (uint32_t i = 0; i < maxVertIndex; ++i) {
			P[i] = verts[i];
			bbox.extendBy(P[i]);
		}

		// allocate memory to store triangle indices
		trisIndex = std::unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);

		// Generate the triangle index array Keep in mind that there is generally 1 vertex attribute for each vertex of each face.
		// So for example if you have 2 quads, you only have 6 vertices but you have 2 * 4 vertex attributes (that is 8 normals,
		// 8 texture coordinates, etc.). So the easiest lazziest method in our triangle mesh, is to create a new array for each
		// supported vertex attribute (st, normals, etc.) whose size is equal to the number of triangles multiplied by 3, and then
		// set the value of the vertex attribute at each vertex of each triangle using the input array (normals[], st[], etc.)
		uint32_t l = 0;
		N = std::unique_ptr<Vec3f []>(new Vec3f[numTris * 3]);
		texCoordinates = std::unique_ptr<Vec2f []>(new Vec2f[numTris * 3]);
		for (uint32_t i = 0, k = 0; i < nfaces; ++i) { // for each face
			for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) { // for each triangle in the face
				trisIndex[l] = vertsIndex[k];
				trisIndex[l + 1] = vertsIndex[k + j + 1];
				trisIndex[l + 2] = vertsIndex[k + j + 2];
				N[l] = normals[k];
				N[l + 1] = normals[k + j + 1];
				N[l + 2] = normals[k + j + 2];
				texCoordinates[l] = st[k];
				texCoordinates[l + 1] = st[k + j + 1];
				texCoordinates[l + 2] = st[k + j + 2];
				l += 3;
			}
			k += faceIndex[i];
		}

		// you can use move if the input geometry is already triangulated
		// N = std::move(normals); // transfer ownership
		// texCoordinates = std::move(st); // transfer ownership
		mailbox.resize(numTris, 0); // should all be initialized with 0
	}

	TriangleMesh(
		const Matrix44f &o2w,
		const uint32_t nfaces,
		const std::unique_ptr<uint32_t []> &faceIndex,
		const std::unique_ptr<uint32_t []> &vertsIndex,
		const std::unique_ptr<Vec3f []> &verts,
		std::unique_ptr<Vec3f []> &normals,
		std::unique_ptr<Vec2f []> &st) :
			Object(o2w),
			numTris(0)
	{
		uint32_t k = 0, maxVertIndex = 0;
		// find out how many triangles we need to create for this mesh
		for (uint32_t i = 0; i < nfaces; ++i) {
			numTris += faceIndex[i] - 2;
			for (uint32_t j = 0; j < faceIndex[i]; ++j) {
				if (vertsIndex[k + j] > maxVertIndex)
					maxVertIndex = vertsIndex[k + j];
			}
			k += faceIndex[i];
		}
		maxVertIndex += 1;

		// allocate memory to store the position of the mesh vertices
		P = std::unique_ptr<Vec3f []>(new Vec3f[maxVertIndex]);
		for (uint32_t i = 0; i < maxVertIndex; ++i) {
			objectToWorld.multVecMatrix(verts[i], P[i]);
			bbox.extendBy(verts[i]);
		}
		objectToWorld.multVecMatrix(bbox[0], bbox[0]);
		objectToWorld.multVecMatrix(bbox[1], bbox[1]);

		// allocate memory to store triangle indices
		trisIndex = std::unique_ptr<uint32_t []>(new uint32_t [numTris * 3]);
		N = std::unique_ptr<Vec3f []>(new Vec3f[numTris * 3]);
		texCoordinates = std::unique_ptr<Vec2f []>(new Vec2f[numTris * 3]);
		mailbox.resize(numTris, 0); // should all be initialized with 0

		Matrix44f transformNormals = worldToObject.transpose();
		// generate the triangle index array and set normals and st coordinates
		uint32_t l = 0;
		for (uint32_t i = 0, k = 0; i < nfaces; ++i) { // for each face
			for (uint32_t j = 0; j < faceIndex[i] - 2; ++j) { // for each triangle in the face
				trisIndex[l] = vertsIndex[k];
				trisIndex[l + 1] = vertsIndex[k + j + 1];
				trisIndex[l + 2] = vertsIndex[k + j + 2];
				transformNormals.multDirMatrix(normals[k], N[l]);
				transformNormals.multDirMatrix(normals[k + j + 1], N[l + 1]);
				transformNormals.multDirMatrix(normals[k + j + 2], N[l + 2]);
				N[l].normalize();
				N[l + 1].normalize();
				N[l + 2].normalize();
				texCoordinates[l] = st[k];
				texCoordinates[l + 1] = st[k + j + 1];
				texCoordinates[l + 2] = st[k + j + 2];
				l += 3;
			}
			k += faceIndex[i];
		}
	}

	bool intersect(const Vec3f& rayOrig, const Vec3f& rayDir, float &tNear) const
	{
		// naive approach, loop over all triangles in the mesh and return true if one
		// of the triangles at least is intersected
		uint32_t j = 0;
		bool intersected = false;
		uint32_t triIndex;
		// tNear should be set infinity first time this function is called and it
		// will get eventually smaller as the ray intersects geometry
		for (uint32_t i = 0; i < numTris; ++i) {
			float t, u, v;
			if (scratch::geometry_utils::rayTriangleIntersect(rayOrig, rayDir,
				P[trisIndex[j]],
				P[trisIndex[j + 1]],
				P[trisIndex[j + 2]], t, u, v) && t < tNear)
			{
				tNear = t;
				triIndex = i;
				intersected |= true;
			}
			j += 3;
		}

		return intersected;
	}

	bool intersect(const Vec3f& rayOrig, const Vec3f& rayDir, float &tNear,
			std::atomic<uint32_t>& triangle_test_counter,
			std::atomic<uint32_t>& triangle_intersection_counter) const
	{
		// naive approach, loop over all triangles in the mesh and return true if one
		// of the triangles at least is intersected
		uint32_t j = 0;
		bool intersected = false;
		uint32_t triIndex;
		// tNear should be set infinity first time this function is called and it
		// will get eventually smaller as the ray intersects geometry
		for (uint32_t i = 0; i < numTris; ++i) {
			float t, u, v;
			triangle_test_counter++;
			if (scratch::geometry_utils::rayTriangleIntersect(rayOrig, rayDir,
				P[trisIndex[j]],
				P[trisIndex[j + 1]],
				P[trisIndex[j + 2]], t, u, v) && t < tNear)
			{
				tNear = t;
				triIndex = i;
				intersected |= true;
				triangle_intersection_counter++;
			}
			j += 3;
		}

		return intersected;
	}

	// Test if the ray intersects this triangle mesh with cross products
	bool intersect(const Vec3f &orig, const Vec3f &dir, float &tNear,
			uint32_t &triIndex,
			Vec2f &uv) const
	{
		uint32_t j = 0;
		bool intersected = false;
		for (uint32_t k = 0; k < numTris; ++k) {
			const Vec3f& v0 = P[trisIndex[j]];
			const Vec3f& v1 = P[trisIndex[j + 1]];
			const Vec3f& v2 = P[trisIndex[j + 2]];
			float t, u, v;
			if (scratch::geometry_utils::rayTriangleIntersect(orig, dir,
				v0, v1, v2, t, u, v) && t < tNear) {
				tNear = t;
				uv.x = u;
				uv.y = v;
				triIndex = k;
				intersected |= true;
			}
			j += 3;
		}

		return intersected;
	}

	void getSurfaceProperties(
		const Vec3f &hitPoint,
		const Vec3f &viewDirection,
		const uint32_t &triIndex,
		const Vec2f &uv,
		Vec3f &hitNormal,
		Vec2f &hitTextureCoordinates) const
	{
		if (smoothShading) {
			// vertex normals
			const Vec3f &n0 = N[triIndex * 3];
			const Vec3f &n1 = N[triIndex * 3 + 1];
			const Vec3f &n2 = N[triIndex * 3 + 2];
			// interpolate with barycentric coordinates
			hitNormal = (1 - uv.x - uv.y) * n0 + uv.x * n1 + uv.y * n2;
		}
		else {
			// face normal - no interpolation
			const Vec3f &v0 = P[trisIndex[triIndex * 3]];
			const Vec3f &v1 = P[trisIndex[triIndex * 3 + 1]];
			const Vec3f &v2 = P[trisIndex[triIndex * 3 + 2]];
			hitNormal = (v1 - v0).crossProduct(v2 - v0);
		}

		hitNormal.normalize();

		// texture coordinates
		const Vec2f &st0 = texCoordinates[triIndex * 3];
		const Vec2f &st1 = texCoordinates[triIndex * 3 + 1];
		const Vec2f &st2 = texCoordinates[triIndex * 3 + 2];
		// interpolate with barycentric coordinates
		hitTextureCoordinates = (1 - uv.x - uv.y) * st0 + uv.x * st1 + uv.y * st2;
	}

	Vec3f evalDiffuseColor(const Vec2f &st) const
	{
		float scale = 5.f;
		float pattern = (fmodf(st.x * scale, 1.f) > 0.5f) ^ (fmodf(st.y * scale, 1.f) > 0.5f);
		return mix(Vec3f(0.815f, 0.235f, 0.031f), Vec3f(0.937f, 0.937f, 0.231f), pattern);
	}

	// member variables
	uint32_t numTris; // number of triangles
	std::unique_ptr<Vec3f []> P; // triangles vertex position
	std::unique_ptr<uint32_t []> trisIndex; // vertex index array
	std::unique_ptr<Vec3f []> N; // triangles vertex normals
	std::unique_ptr<Vec2f []> texCoordinates; // triangles texture coordinates

	// [comment]
	// Mailboxes are used by the Grid acceleration structure
	// [/comment]
	mutable std::vector<uint32_t> mailbox;
};

// [comment]
// Compute the position of a point along a Bezier curve at t [0:1]
// [/comment]
Vec3f evalBezierCurve(const Vec3f *P, const float &t)
{
	float b0 = (1 - t) * (1 - t) * (1 - t);
	float b1 = 3 * t * (1 - t) * (1 - t);
	float b2 = 3 * t * t * (1 - t);
	float b3 = t * t * t;

	return P[0] * b0 + P[1] * b1 + P[2] * b2 + P[3] * b3;
}

Vec3f evalBezierPatch(const Vec3f *controlPoints, const float &u, const float &v)
{
	Vec3f uCurve[4];
	for (int i = 0; i < 4; ++i) {
		uCurve[i] = evalBezierCurve(controlPoints + 4 * i, u);
	}

	return evalBezierCurve(uCurve, v);
}

Vec3f derivBezier(const Vec3f *P, const float &t)
{
	return -3 * (1 - t) * (1 - t) * P[0] +
		(3 * (1 - t) * (1 - t) - 6 * t * (1 - t)) * P[1] +
		(6 * t * (1 - t) - 3 * t * t) * P[2] +
		3 * t * t * P[3];
}

// [comment]
// Compute the derivative of a point on Bezier patch along the u parametric direction
// [/comment]
Vec3f dUBezier(const Vec3f *controlPoints, const float &u, const float &v)
{
	Vec3f P[4];
	Vec3f vCurve[4];
	for (int i = 0; i < 4; ++i) {
		P[0] = controlPoints[i];
		P[1] = controlPoints[4 + i];
		P[2] = controlPoints[8 + i];
		P[3] = controlPoints[12 + i];
		vCurve[i] = evalBezierCurve(P, v);
	}

	return derivBezier(vCurve, u);
}

// [comment]
// Compute the derivative of a point on Bezier patch along the v parametric direction
// [/comment]
Vec3f dVBezier(const Vec3f *controlPoints, const float &u, const float &v)
{
	Vec3f uCurve[4];
	for (int i = 0; i < 4; ++i) {
		uCurve[i] = evalBezierCurve(controlPoints + 4 * i, u);
	}

	return derivBezier(uCurve, v);
}
