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

#include "geometry.h"

enum RayType {
	kPrimaryRay,
	kShadowRay
};

class LightLite
{
public:
	LightLite(const Vec3f &p, const Vec3f &i) : position(p), intensity(i) {}
	Vec3f position;
	Vec3f intensity;
};

class Light
{
public:
	Light(const Matrix44f &l2w, const Vec3f &c = 1.f, const float &i = 1.f) : lightToWorld(l2w), color(c), intensity(i) {}
	Light(const Vec3f& c, const float& i) : color(c), intensity(i) {}
	virtual ~Light() {}
	virtual void illuminate(const Vec3f &P, Vec3f &, Vec3f &, float &) const = 0;
	Vec3f color;
	float intensity;
	Matrix44f lightToWorld;
};

class DistantLight : public Light
{
public:
	DistantLight(const Matrix44f &l2w, const Vec3f &c = 1.f, const float &i = 1.f) : Light(l2w, c, i)
	{
		l2w.multDirMatrix(Vec3f(0.f, 0.f, -1.f), dir);
		dir.normalize(); // in case the matrix scales the light
	}

	void illuminate(const Vec3f &P, Vec3f &lightDir, Vec3f &lightIntensity, float &distance) const
	{
		lightDir = dir;
		lightIntensity = color * intensity;
		distance = kInfinity;
	}

	Vec3f dir;
};

class PointLight : public Light
{
public:
	PointLight(const Matrix44f &l2w, const Vec3f &c = 1.f, const float &i = 1.f) : Light(l2w, c, i)
	{
		l2w.multVecMatrix(Vec3f(0), pos);
	}

	PointLight(const Vec3f& p, const Vec3f& c, const float& i) : Light(c, i), pos(p) {}

	// P: is the shaded point
	void illuminate(const Vec3f &P, Vec3f &lightDir, Vec3f &lightIntensity, float &distance) const
	{
		lightDir = (P - pos);
		float r2 = lightDir.norm();
		distance = sqrt(r2);
		lightDir.x /= distance, lightDir.y /= distance, lightDir.z /= distance;
		// avoid division by 0
		lightIntensity = color * intensity / (4 * M_PI * r2);
	}

	Vec3f pos;
};
