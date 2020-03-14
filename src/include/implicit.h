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

// [comment]
// ImplicitShape base class
// [/comment]
class ImplicitShape
{
public:
	virtual float getDistance(const Vec3f& from) const = 0;
	virtual ~ImplicitShape() {}
};

// [comment]
// Implicit sphere surface
// [/comment]
class ImplicitSphere : public ImplicitShape
{
public:
	ImplicitSphere(const Vec3f& c, const float& r) : center(c), radius(r) {}
	float getDistance(const Vec3f& from) const
	{
		return (from - center).length() - radius;
	}
	float radius;
	Vec3f center;
};

// [comment]
// Implicit plane surface
// [/comment]
class ImplicitPlane : public ImplicitShape
{
public:
	ImplicitPlane(const Vec3f& nn = Vec3f(0, 1, 0), const Vec3f& pp = Vec3f(0)) :
		n(nn), pointOnPlane(pp) {}
	float getDistance(const Vec3f& from) const
	{
		return n.dotProduct(from - pointOnPlane);
	}
	Vec3f n, pointOnPlane;
};

// [comment]
// Implicit torus surface
// [/comment]
class ImplicitTorus : public ImplicitShape
{
public:
	ImplicitTorus(const float& _r0, const float& _r1) : r0(_r0), r1(_r1) {}
	float getDistance(const Vec3f& from) const
	{
		// reduce 3D point to 2D point
		float tmpx = sqrtf(from.x * from.x + from.z * from.z) - r0;
		float tmpy = from.y;

		// distance to cicle
		return sqrtf(tmpx * tmpx + tmpy * tmpy) - r1;
	}
	float r0, r1;
};

// [comment]
// Implicit cube surface
// [/comment]
class ImplicitCube : public ImplicitShape
{
public:
	ImplicitCube(const Vec3f &_corner) : corner(_corner) {}
	float getDistance(const Vec3f& from) const
	{
#if 0
		// first transform the input point into the object's "object-space".
		float scale = 2.f;

		// this matrix doesn't scale the object
		Matrix44f objectToWorld(0.542903, -0.545887, 0.638172, 0, 0.778733, 0.611711, -0.139228, 0, -0.314374, 0.572553, 0.7572, 0, 0, 1.459974, 0, 1);
		Matrix44f worldToObject = objectToWorld.inverse();

		Vec3f fromObjectSpace = from;
		worldToObject.multVecMatrix(from, fromObjectSpace);
#else
		Vec3f fromObjectSpace = from;
		float scale = 1;
#endif
		fromObjectSpace *= 1.f / scale;
		fromObjectSpace.x = std::fabs(fromObjectSpace.x);
		fromObjectSpace.y = std::fabs(fromObjectSpace.y);
		fromObjectSpace.z = std::fabs(fromObjectSpace.z);

		// now compute the distance from the point to the neares point on the surface of the object
		Vec3f d = fromObjectSpace - corner;

		Vec3f dmax = d;
		dmax.x = std::max(dmax.x, 0.f);
		dmax.y = std::max(dmax.y, 0.f);
		dmax.z = std::max(dmax.z, 0.f);

		// don't forget to apply the scale back
		return scale * (std::min(std::max(d.x, std::max(d.y, d.z)), 0.f) + dmax.length());
	}
	Vec3f corner{0.25, 0.25, 0.25};
};

struct unionFunc
{
	float operator() (float a, float b) const { return std::min(a, b); }
};

struct subtractFunc
{
	float operator() (float a, float b) const { return std::max(-a, b); }
};

struct intersectionFunc
{
	float operator() (float a, float b) const { return std::max(a, b); }
};

struct blendFunc
{
	blendFunc(const float &_k) : k(_k) {}
	float operator() (float a, float b) const
	{
		float res = exp(-k * a) + exp(-k * b);
		return -log(std::max(0.0001f, res)) / k;
	}
	float k;
};

struct mixFunc
{
	mixFunc(const float &_t) : t(_t) {}
	float operator() (float a, float b) const
	{
		return a * (1 -t) + b * t;
	}
	float t;
};

// [comment]
// Combining implict shapes using CSG
// [/comment]
template<typename Op, typename ... Args>
class CSG : public ImplicitShape
{
public:
	CSG(
		const std::shared_ptr<ImplicitShape> s1,
		const std::shared_ptr<ImplicitShape> s2,
		Args&& ... args) : op(std::forward<Args>(args) ...), shape1(s1), shape2(s2)
	{}
	float getDistance(const Vec3f& from) const
	{
		return op(shape1->getDistance(from), shape2->getDistance(from));
	}
	Op op;
	const std::shared_ptr<ImplicitShape> shape1, shape2;
};

// [comment]
// Blobbies
// [/comment]
class SoftObject : public ImplicitShape
{
	struct Blob
	{
		float R; // radius
		Vec3f c; // blob center
	};
public:
	SoftObject()
	{
	#if 1
		blobbies.push_back({2.0, Vec3f(-1, 0, 0)});
		blobbies.push_back({1.5, Vec3f( 1, 0, 0)});
	#else
		for (size_t i = 0; i < 20; ++i) {
			float radius = (0.3 + drand48() * 1.3);
			Vec3f c((0.5 - drand48()) * 3, (0.5 - drand48()) * 3, (0.5 - drand48()) * 3);
			blobbies.push_back({radius, c});
		}
	#endif
	}
	float getDistance(const Vec3f& from) const
	{
		float sumDensity = 0;
		float sumRi = 0;
		float minDistance = kInfinity;
		for (const auto& blob: blobbies) {
			float r = (blob.c - from).length();
			if (r <=  blob.R) {
				// this can be factored for speed if you want
				sumDensity += 2 * (r * r * r) / (blob.R * blob.R * blob.R) -
    			3 * (r * r) / (blob.R * blob.R) + 1;
			}
			minDistance = std::min(minDistance, r - blob.R);
			sumRi += blob.R;
		}

		return std::max(minDistance, (magic - sumDensity) / ( 3 / 2.f * sumRi));
	}
	float magic{ 0.2 };
	std::vector<Blob> blobbies;
};

using Union = CSG<unionFunc>;
using Subtract = CSG<subtractFunc>;
using Intersect = CSG<intersectionFunc>;
using Blend = CSG<blendFunc, float>;
using Mix = CSG<mixFunc, float>;