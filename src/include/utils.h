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

#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI (3.14159265358979323846f)
#endif

constexpr float kInfinity = std::numeric_limits<float>::max();
constexpr float kEpsilon = 1e-8;

namespace scratch
{
namespace utils
{
	inline float mix(const float &a, const float &b, const float &mix)
	{
		return b * mix + a * (1 - mix);
	}

	inline float min3(const float &a, const float &b, const float &c)
	{
		return std::min(a, std::min(b, c));
	}

	inline float max3(const float &a, const float &b, const float &c)
	{
		return std::max(a, std::max(b, c));
	}

	inline float clamp(const float &lo, const float &hi, const float &v)
	{
		return std::max(lo, std::min(hi, v));
	}

	template<typename T = float>
	inline T clamp(const T &v, const T &lo, const T &hi)
	{
		return std::max(lo, std::min(hi, v));
	}

	inline float deg2rad(const float &deg)
	{
		return deg * M_PI / 180;
	}

	inline float modulo(const float &f)
	{
		return f - std::floor(f);
	}

	template<typename T = float>
	inline float edgeFunction(const T (&a)[2], const T (&b)[2], const T (&c)[2])
	{
		return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
	}

	template<typename T = float>
	inline float edgeFunction(const T (&a)[3], const T (&b)[3], const T (&c)[3])
	{
		return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
	}

	template<typename T = float>
	inline T lerp(const T &lo, const T &hi, const T &t)
	{
		return lo * (1 - t) + hi * t;
	}

	template<typename T = float>
	inline T smoothstep(const T &t)
	{
		return t * t * (3 - 2 * t);
	}

	inline bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
	{
		float discr = b * b - 4 * a * c;
		if (discr < 0) {
			return false;
		} else if (discr == 0) {
			x0 = x1 = - 0.5f * b / a;
		} else {
			float q = (b > 0) ?
			-0.5f * (b + sqrt(discr)) :
			-0.5f * (b - sqrt(discr));
			x0 = q / a;
			x1 = c / q;
		}
		if (x0 > x1) {
			std::swap(x0, x1);
		}
		return true;
	}
}
}